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
  * The set is an associated array mapping an event_id (an event type) to a throttler instance.
  */
template <typename T>
class JfrEventThrottlerSet : public JfrCHeapObj {
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

  ~JfrEventThrottlerSet() {
    for (int i = 0; i <= LAST_EVENT_ID; i++) {
      delete _set[i];
      _set[i] = NULL;
    }
  }

  T* map(JfrEventId event_id) {
    return _set[event_id];
  }
};

static JfrEventThrottlerSet<JfrEventThrottler>* _throttlers = NULL;
static const uint64_t default_window_duration_ms = 200;
static const uint64_t default_amortization_count = 1; // in number of windows
static const JfrSamplerParams _disabled_params = { 0, 0, 0 };
// The (low) rate, for which the overhead of context switching windows
// with durations lower than one seconds is not justified.
static const uint64_t max_rate_for_per_second_window = 9;

JfrEventThrottler::JfrEventThrottler(JfrEventId event_id) :
  JfrAdaptiveSampler(),
  _event_id(event_id),
  _last_params(_disabled_params),
  _last_rate_limit_per_second(0) {}

bool JfrEventThrottler::initialize() {
  return JfrAdaptiveSampler::initialize();
}

bool JfrEventThrottler::create() {
  assert(_throttlers == NULL, "invariant");
  _throttlers = new JfrEventThrottlerSet<JfrEventThrottler>();
  return _throttlers != NULL && _throttlers->initialize();
}

void JfrEventThrottler::destroy() {
  delete _throttlers;
  _throttlers = NULL;
}

inline void samples_per_window(JfrSamplerParams& params, uint64_t rate_limit_per_second, const JfrSamplerWindow* expired) {
  assert(rate_limit_per_second > 0, "invariant");
  assert(expired != NULL, "invariant");
  if (rate_limit_per_second <= max_rate_for_per_second_window) {
    params.samples_per_window = rate_limit_per_second;
    return;
  }
  // window duration is in milliseconds and rate_limit in samples per second
  const double rate_limit_per_ms = static_cast<double>(rate_limit_per_second) / static_cast<double>(MILLIUNITS);
  params.samples_per_window = floor(rate_limit_per_ms * default_window_duration_ms);
}

inline void window_duration(JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.samples_per_window > 0, "invariant");
  if (params.samples_per_window <= max_rate_for_per_second_window) {
    params.window_duration_ms = MILLIUNITS; // 1 second
    return;
  }
  params.window_duration_ms = default_window_duration_ms;
}

inline void amortization(JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  params.amortization_windows = default_amortization_count;
}

inline void update_params(JfrSamplerParams& params, uint64_t rate_limit_per_second, const JfrSamplerWindow* expired) {
  samples_per_window(params, rate_limit_per_second, expired);
  window_duration(params, expired);
  amortization(params, expired);
}

/*
 * This function is the control mechansism for the AdaptiveSampler engine.
 * The engine calls this function when a sampler window has expired, providing the
 * client with an opportunity to perform some analysis. To reciprocate the favor,
 * the client returns an updated set of parameters for the engine to apply to
 * the next window. Try to keep this relatively quick, since the engine
 * is currently inside a critical section, in the process of context switching windows.
 */
JfrSamplerParams JfrEventThrottler::next_window_params(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  const int64_t rate_limit_per_second = JfrEventSetting::ratelimit(_event_id);
  if (rate_limit_per_second == 0) {
    return _disabled_params;
  }
  if (rate_limit_per_second == _last_rate_limit_per_second) {
    return _last_params;
  }
  // update_params modifies _last_params in-place
  update_params(_last_params, rate_limit_per_second, expired);
  _last_rate_limit_per_second = rate_limit_per_second;
  return _last_params;
}

bool JfrEventThrottler::sample(int64_t timestamp /* 0 */) {
  return JfrEventSetting::ratelimit(_event_id) != 0 ? JfrAdaptiveSampler::sample(timestamp) : true;
}

bool JfrEventThrottler::sample(JfrEventId event_id, int64_t timestamp) {
  assert(_throttlers != NULL, "JfrEventThrottler has not been properly initialized");
  JfrEventThrottler* const throttler = _throttlers->map(event_id);
  return throttler != NULL ? throttler->sample(timestamp) : true;
}
