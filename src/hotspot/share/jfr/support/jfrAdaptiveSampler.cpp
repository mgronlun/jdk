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
#include "jfr/support/jfrAdaptiveSampler.hpp"
#include "jfr/utilities/jfrTime.hpp"
#include "jfr/utilities/jfrTimeConverter.hpp"
#include "jfr/utilities/jfrTryLock.hpp"
#include "runtime/atomic.hpp"
#include "runtime/interfaceSupport.inline.hpp"
#include "runtime/samplerSupport.hpp"
#include "utilities/globalDefinitions.hpp"
#include <cmath>

inline int64_t now() {
  return JfrTicks::now().value();
}

inline int64_t millis_to_countertime(int64_t millis) {
  return JfrTimeConverter::nanos_to_countertime(millis * 1000000);
}

class AdaptiveSampler::Window : public JfrCHeapObj {
  friend class AdaptiveSampler;
 private:
  SamplerWindowParams _params = { -1, -1 };
  int64_t _start_ticks;
  int64_t _end_ticks;
  double _probability;
  size_t _samples_budget;
  volatile size_t _running_count;
  volatile size_t _sample_count;
  Window(double probablity, size_t samples_budget);
  void reinitialize(double probability, size_t samples_budget, SamplerWindowParams params);
  bool sample_attempt();
  bool random_trial();

  size_t sample_count() const {
    const size_t count = Atomic::load(&_sample_count);
    return count < _samples_budget ? count : _samples_budget;
  }

  size_t total_count() const {
    return Atomic::load(&_running_count);
  }

 public:
  Window(SamplerWindowParams params, double probability, size_t samples_budget);
  bool should_sample(int64_t timestamp, bool* is_expired);


  const SamplerWindowParams params() const {
    return _params;
  }

  /**
   * Ratio between the requested and the measured window duration
   */
   double adjustment_factor() const {
    return (double)(_end_ticks - _start_ticks) / (now() - _start_ticks);
  }
  
  double adjustment_factor(jlong window_duration_ms) const {
    return (double)millis_to_countertime(window_duration_ms) / (now() - _start_ticks);
  }

  bool is_expired() const {
    return is_expired(now());
  }

  bool is_expired(int64_t timestamp) const {
    return timestamp >= _end_ticks;
  }

  double probability() const {
    return Atomic::load_acquire(&_probability);
  }
};

AdaptiveSampler::Window::Window(SamplerWindowParams params, double probability, size_t samples_budget) :
    _probability(probability),
    _params(params),
    _samples_budget(samples_budget),
    _start_ticks(params.duration == -1 ? 0 : now()),
    _end_ticks(params.duration == -1 ? 0 : _start_ticks + millis_to_countertime(params.duration)),
    _running_count(0),
    _sample_count(0) {}

AdaptiveSampler::Window::Window(double probability, size_t samples_budget) :
    _probability(probability),
    _samples_budget(samples_budget),
    _start_ticks(_params.duration == -1 ? 0 : now()),
    _end_ticks(_params.duration == -1 ? 0 : _start_ticks + millis_to_countertime(_params.duration)),
    _running_count(0),
    _sample_count(0) {}

inline void AdaptiveSampler::Window::reinitialize(double probability, size_t samples_budget, SamplerWindowParams params) {
  _probability = probability;
  Atomic::store(&_running_count, (size_t)0);
  Atomic::store(&_sample_count, (size_t)0);
  _params = params;
  _samples_budget = samples_budget;
  _start_ticks = params.duration == -1 ? 0 : now(),
  _end_ticks = params.duration == -1 ? 0 : _start_ticks + millis_to_countertime(params.duration);
}

static SamplerSupport* sampler_support() {
  static THREAD_LOCAL SamplerSupport* _sampler_support = NULL;
  if (_sampler_support == NULL) {
    _sampler_support = new SamplerSupport(true);
  }
  return _sampler_support;
}

inline double next_random_uniform() {
  return sampler_support()->next_random_uniform();
}

/*
 * The adaptive sampler equips each window with a target ratio,
 * which describes the relation between the number of samples to be taken
 * over the estimated number of sample attempts. This ratio, or fraction,
 * is expressed as a probability value.
 *
 * The should_sample() predicate can therefore be modelled as a simple random trial,
 * the outcome of which decides if a sample should be taken or not.
 */

// P(X = x) refers to the probability that the random variable X
// is less than a particular value x.
bool AdaptiveSampler::Window::random_trial() {
  const double x = probability();
  if (x < 1) {
    // let's see what the random variable X is
    const double X = next_random_uniform();
    if (X >= x) {
      // Event did not occur.
      return false;
    }
  }
  // Event did occur. Select for sampling if still within the samples budget.
  return Atomic::add(&_sample_count, (size_t)1, memory_order_acq_rel) <= _samples_budget;
}

bool AdaptiveSampler::Window::sample_attempt() {
  // Increment number of attempts
  Atomic::add(&_running_count, (size_t)1, memory_order_acq_rel);
  return random_trial();
}

bool AdaptiveSampler::Window::should_sample(int64_t timestamp, bool* expired_window) {
  assert(expired_window != NULL, "invariant");
  *expired_window = is_expired(timestamp);
  return *expired_window ? false : sample_attempt();
}

inline double compute_interval_alpha(size_t interval) {
  return 1 - std::pow(interval, (double)-1 / (double)interval);
}

AdaptiveSampler::AdaptiveSampler(size_t window_lookback_cnt, size_t budget_lookback_cnt) :
    _window_0(NULL),
    _window_1(NULL),
    _active_window(NULL),
    _window_lookback_alpha(compute_interval_alpha(window_lookback_cnt)),
    _budget_lookback_alpha(compute_interval_alpha(budget_lookback_cnt)),
    _samples_budget(-1 * (1 + budget_lookback_cnt)),
    _probability(1),
    _budget_lookback_cnt(budget_lookback_cnt),
    _avg_samples(std::numeric_limits<double>::quiet_NaN()),
    _avg_count(0),
    _lock(0) {}

bool AdaptiveSampler::initialize() {
  assert(_window_0 == NULL, "invariant");
  _window_0 = new Window(_probability, (size_t)_samples_budget);
  if (_window_0 == NULL) {
    return false;
  }
  assert(_window_1 == NULL, "invariant");
  _window_1 = new Window(_probability, (size_t)_samples_budget);
  if (_window_1 == NULL) {
    return false;
  }
  _active_window = _window_0;
  return true;
}

AdaptiveSampler::~AdaptiveSampler() {
  delete _window_0;
  delete _window_1;
}

bool AdaptiveSampler::should_sample() {
  bool expired_window = false;
  int64_t timestamp = now();
  bool sample = active_window()->should_sample(timestamp, &expired_window);
  if (expired_window) {
    {
      JfrTryLock rotation_lock(&_lock);
      if (rotation_lock.acquired()) {
        rotate_window();
      }
    }
    timestamp = now();
    sample = active_window()->should_sample(timestamp, &expired_window);
  }
  return sample;
}

inline AdaptiveSampler::Window* AdaptiveSampler::active_window() {
  return Atomic::load_acquire(&_active_window);
}

// Called when the previous window has expired
// and exlusively by the holder of the try_lock
void AdaptiveSampler::rotate_window() {
  assert(_lock, "invariant");
  Window* const prev_window = active_window();
  if (!prev_window->is_expired()) {
    // someone took care of it
    return;
  }
  SamplerWindowParams new_params = new_window_params();
  recalculate_averages(prev_window, new_params);
  recalculate_probability(new_params);
  install_new_window(prev_window, new_params);
  assert(prev_window->is_expired(), "invariant");
}

void AdaptiveSampler::recalculate_averages(const AdaptiveSampler::Window* prev_window, SamplerWindowParams new_params) {
  assert(prev_window != NULL, "invariant");
  assert(prev_window == active_window(), "invariant");
  const SamplerWindowParams previous_params = prev_window->params();
  bool is_dummy = previous_params.duration == -1;
  const double adjustment_factor = is_dummy ? prev_window->adjustment_factor(new_params.duration) : prev_window->adjustment_factor();
  const double samples = prev_window->sample_count() * adjustment_factor;
  const double total_count = prev_window->total_count() * adjustment_factor;

  if (!is_dummy) {
    _avg_samples = std::isnan(_avg_samples) ? samples : _avg_samples + _budget_lookback_alpha * (samples - _avg_samples);
  }
  _samples_budget = fmax<double>(new_params.sample_count - _avg_samples, 0) * _budget_lookback_cnt;

  // fprintf(stdout, "=== avg samples: %f, samples: %f, adjustment: %f\n", _avg_samples, samples, adjustment_factor);

  if (_avg_count == 0) {
    _avg_count = total_count;
  } else {
    // need to convert int '*_count' variables to double to prevent bit overflow
    _avg_count = _avg_count + _window_lookback_alpha * ((double)total_count - (double)_avg_count);
  }
}

void AdaptiveSampler::recalculate_probability(SamplerWindowParams params) {
  if (_avg_count == 0) {
    _probability = 1;
    return;
  }
  const double p = (params.sample_count + _samples_budget) / (double)_avg_count;
  _probability = p < 1 ? p : 1;
}

void AdaptiveSampler::install_new_window(AdaptiveSampler::Window* prev_window, SamplerWindowParams new_params) {
  assert(prev_window != NULL, "invariant");
  assert(prev_window == active_window(), "invariant");
  Window* const next_window = idle_window(prev_window);
  assert(next_window != prev_window, "invariant");
  next_window->reinitialize(_probability, _samples_budget, new_params);
  Atomic::release_store(&_active_window, next_window);
}

inline AdaptiveSampler::Window* AdaptiveSampler::idle_window(const AdaptiveSampler::Window* prev_window) {
  assert(prev_window != NULL, "invariant");
  return prev_window == _window_0 ? _window_1 : _window_0;
}

