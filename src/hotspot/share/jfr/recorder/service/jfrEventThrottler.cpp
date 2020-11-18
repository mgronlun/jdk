/*
 * Copyright (c) 2020, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 *
 */

#include "precompiled.hpp"
#include "jfr/recorder/jfrEventSetting.inline.hpp"
#include "jfr/recorder/service/jfrEventThrottler.hpp"
#include "jfr/utilities/jfrAllocation.hpp"
#include "runtime/semaphore.inline.hpp"

 /*
  * This is an associative array.
  * It maps an event type to its throttler instance using the event id as the key.
  */
template <typename T>
class JfrEventThrottlerMap : public JfrCHeapObj {
 private:
  T* _set[LAST_EVENT_ID + 1];
 public:
  bool initialize() {
    for (int i = FIRST_EVENT_ID; i <= LAST_EVENT_ID; i++) {
      _set[i] = new T(static_cast<JfrEventId>(i));
      if (_set[i] == NULL || !_set[i]->initialize()) {
        return false;
      }
    }
    return true;
  }

  ~JfrEventThrottlerMap() {
    for (int i = FIRST_EVENT_ID; i <= LAST_EVENT_ID; i++) {
      delete _set[i];
      _set[i] = NULL;
    }
  }

  T* get(JfrEventId event_id) {
    return _set[event_id];
  }
};

class JfrEventThrottlerSettingsLocker : public StackObj {
 private:
  Semaphore* const _semaphore;
 public:
  JfrEventThrottlerSettingsLocker(JfrEventThrottler* throttler) : _semaphore(throttler->_semaphore) {
    assert(_semaphore != NULL, "invariant");
    _semaphore->wait();
  }
  ~JfrEventThrottlerSettingsLocker() {
    _semaphore->signal();
  }
};

static const JfrSamplerParams _disabled_params = {
                                                   JfrSamplerParams::unused,
                                                   JfrSamplerParams::unused,
                                                   JfrSamplerParams::unused,
                                                   0,
                                                   0,
                                                   false
                                                 };

JfrEventThrottler::JfrEventThrottler(JfrEventId event_id) :
  JfrAdaptiveSampler(),
  _semaphore(NULL),
  _event_id(event_id),
  _last_params(_disabled_params),
  _disabled(false),
  _update(false) {}

JfrEventThrottler::~JfrEventThrottler() {
  delete _semaphore;
}

bool JfrEventThrottler::initialize() {
  _semaphore = new Semaphore(1);
  return _semaphore != NULL && JfrAdaptiveSampler::initialize();
}

/*
 * The event throttler can be configured in three dimensions (all optional):
 *
 * - rate per second  throttle dynamically to maintain a continuous, maximal rate per second
 * - probability      throttle using a probability
 * - nth selection    throttle by selecting every nth sample point
 *
 * If a rate is specified, the probability and / or nth selection becomes relative to the rate.
 */
bool JfrEventThrottler::configure(intptr_t rate_per_second, double probability, intptr_t nth_selection) {
  JfrEventThrottlerSettingsLocker sl(this);
  _rate_per_second = rate_per_second;
  _probability = probability;
  _nth_selection = nth_selection;
  _update = true;
  return true;
}

static JfrEventThrottlerMap<JfrEventThrottler>* _throttlers = NULL;

bool JfrEventThrottler::create() {
  assert(_throttlers == NULL, "invariant");
  _throttlers = new JfrEventThrottlerMap<JfrEventThrottler>();
  return _throttlers != NULL && _throttlers->initialize();
}

void JfrEventThrottler::destroy() {
  delete _throttlers;
  _throttlers = NULL;
}

JfrEventThrottler* JfrEventThrottler::for_event(JfrEventId event_id) {
  assert(_throttlers != NULL, "JfrEventThrottler has not been properly initialized");
  return _throttlers->get(event_id);
}

bool JfrEventThrottler::accept(JfrEventId event_id, int64_t timestamp) {
  JfrEventThrottler* const throttler = for_event(event_id);
  assert(throttler != NULL, "invariant");
  return throttler->sample(timestamp);
  //return throttler->_disabled ? true : throttler->sample(timestamp);
}

const intptr_t event_throttler_disabled = -2;
static const size_t default_window_duration_ms = 200;   // 5 windows per second

/*
 * Rates lower than or equal to the 'low rate upper bound', are considered special.
 * They will use a window of duration one second, because the rates are so low they
 * do not justify the overhead of more frequent window rotations.
 */
static const size_t low_rate_upper_bound = 9;

/*
 * Breaks down an overall rate per second to a number of sample points per window.
 */
inline void set_sample_points_per_window(JfrSamplerParams& params, int64_t rate_per_second, const JfrSamplerWindow* expired) {
  assert(rate_per_second != event_throttler_disabled, "invariant");
  assert(expired != NULL, "invariant");
  if (rate_per_second <= low_rate_upper_bound) {
    params.sample_points_per_window = rate_per_second;
    return;
  }
  // Window duration is in milliseconds and the rate_is in sample points per second.
  const double rate_per_ms = static_cast<double>(rate_per_second) / static_cast<double>(MILLIUNITS);
  params.sample_points_per_window = floor(rate_per_ms * default_window_duration_ms);
}

/*
 * The window_lookback_count states the history in number of windows to take into account
 * when the JfrAdaptiveSampler engine is calcualting an expected weigthed moving average (EWMA).
 * It only applies to contexts where a rate is specified. Technically, it determines the alpha
 * coefficient in an EMWA formula.
 */
static const size_t default_window_lookback_count = 25; // 25 windows == 5 seconds (for default window duration)

inline void set_window_duration(JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  if (params.sample_points_per_window != JfrSamplerParams::unused && params.sample_points_per_window > low_rate_upper_bound) {
    params.window_duration_ms = default_window_duration_ms;
    params.window_lookback_count = default_window_lookback_count; // 5 seconds
    return;
  }
  // Lower rates.
  params.window_duration_ms = MILLIUNITS; // 1 second
  params.window_lookback_count = 5; // 5 windows == 5 seconds
}

const JfrSamplerParams& JfrEventThrottler::last_params(const JfrSamplerWindow* expired) {
  JfrEventThrottlerSettingsLocker sl(this);
  if (event_throttler_disabled == _rate_per_second) {
    // The event throttler is "off".
    _disabled = true;
    return _disabled_params;
  }
  _last_params.probability = _probability;
  _last_params.nth_selection = _nth_selection;
  set_sample_points_per_window(_last_params, _rate_per_second, expired);
  set_window_duration(_last_params, expired);
  _last_params.reconfigure = true;
  _update = false;
  return _last_params;
}

/*
 * This is the feedback control loop when using the JfrAdaptiveSampler engine.
 *
 * The engine calls this when a sampler window has expired, providing the
 * client with an opportunity to perform some analysis. To reciprocate, the client
 * returns a set of parameters, possibly updated, for the engine to apply to the next window.
 *
 * Try to keep relatively quick, since the engine is currently inside a critical section,
 * in the process of rotating windows.
 */
const JfrSamplerParams& JfrEventThrottler::next_window_params(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  if (_update) {
    return last_params(expired); // Updates _last_params in-place.
  }
  return _disabled ? _disabled_params : _last_params;
}
