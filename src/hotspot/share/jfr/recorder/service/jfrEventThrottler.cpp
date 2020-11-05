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

static const JfrSamplerParams _disabled_params = { 0, 0, 0 };
static const size_t default_window_duration_ms = 200;
static const size_t window_lookback_count = 25;

JfrEventThrottler::JfrEventThrottler(JfrEventId event_id) :
  JfrAdaptiveSampler(window_lookback_count),
  _event_id(event_id),
  _last_params(_disabled_params),
  _last_rate_limit_per_second(0) {}

bool JfrEventThrottler::initialize() {
  return JfrAdaptiveSampler::initialize();
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
  if (JfrEventSetting::ratelimit(event_id) == 0) {
    return true;
  }
  JfrEventThrottler* const throttler = for_event(event_id);
  return throttler != NULL ? throttler->sample(timestamp) : true;
}

/*
 * Rates lower than or equal to the 'low rate upper bound', are considered special.
 * They will use a window with duration one second, because the rates are so low they
 * do not justify the overhead of more frequent context switching of windows.
 */
static const size_t low_rate_upper_bound = 9;

inline void samples_per_window(JfrSamplerParams& params, uint64_t rate_limit_per_second, const JfrSamplerWindow* expired) {
  assert(rate_limit_per_second > 0, "invariant");
  assert(expired != NULL, "invariant");
  if (rate_limit_per_second <= low_rate_upper_bound) {
    params.sample_points_per_window = rate_limit_per_second;
    return;
  }
  // Window duration is in milliseconds and the rate_limit is in samples per second.
  const double rate_limit_per_ms = static_cast<double>(rate_limit_per_second) / static_cast<double>(MILLIUNITS);
  params.sample_points_per_window = floor(rate_limit_per_ms * default_window_duration_ms);
}

inline void window_duration(JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window > 0, "invariant");
  if (params.sample_points_per_window <= low_rate_upper_bound) {
    params.window_duration_ms = MILLIUNITS; // 1 second
    return;
  }
  params.window_duration_ms = default_window_duration_ms;
}

/*
 * This parameter controls how fast 'accumulated debt' should be amortized.
 * Debt can be thought of as a cumulative error term, and is indicative for how much the sampler
 * is deviating from a set point. Due to the adaptive nature of the sampler, debt accumulates naturally
 * over time, as a function of under-and/or oversampled windows, because of system fluctuations.
 * This parameter controls how fast to catch up / slow down, with the rate given in number of windows.
 * The default is one window, meaning the sampler should attempt to amortize an accumulated debt immediately
 * in the next window, if possible. A low amoritization count will on average keep the sampler
 * closer to the set point. A potential downside with a low amortization count is that some randomness
 * might be sacrificed in the process of having the sampler catch up more quickly.
 */
static const uint64_t default_amortization_count = 1;

inline void amortization(JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  params.amortization_window_count = default_amortization_count;
}

const JfrSamplerParams& JfrEventThrottler::last_params(int64_t rate_limit_per_second, const JfrSamplerWindow* expired) {
  samples_per_window(_last_params, rate_limit_per_second, expired);
  window_duration(_last_params, expired);
  amortization(_last_params, expired);
  _last_rate_limit_per_second = rate_limit_per_second;
  return _last_params;
}

/*
 * This function is the control mechanism for the JfrAdaptiveSampler engine.
 * The engine calls this function when a sampler window has expired, providing the
 * client with an opportunity to perform some analysis. To reciprocate, the client
 * returns an updated set of parameters for the engine to apply to the next window.
 *
 * Try to keep relatively quick, since the engine is currently inside a critical section,
 * in the process of context switching windows.
 */
const JfrSamplerParams& JfrEventThrottler::next_window_params(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  const int64_t rate_limit_per_second = JfrEventSetting::ratelimit(_event_id);
  if (rate_limit_per_second == 0) {
    return _disabled_params;
  }
  if (rate_limit_per_second == _last_rate_limit_per_second) {
    return _last_params;
  }
  // last_params() modifies _last_params in-place
  return last_params(rate_limit_per_second, expired);
}
