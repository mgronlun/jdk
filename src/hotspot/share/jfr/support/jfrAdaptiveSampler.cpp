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

/*
inline
const size_t total_count() {
  return Atomic::load(&_running_count);
}
*/

static SamplerSupport* _sampler_support = NULL;

static SamplerSupport* get_sampler_support() {
  assert(_sampler_support != NULL, "invariant");
  return _sampler_support;
}

inline int64_t millis_to_countertime(int64_t millis) {
  return JfrTimeConverter::nanos_to_countertime(millis * 1000000);
}

inline int64_t now() {
  return JfrTicks::now().value();
}

class AdaptiveSampler::Window : public JfrCHeapObj {
  friend class AdaptiveSampler;
 private:
  const SamplerWindowParams _params = { -1, -1 };
  const int64_t _start_ticks;
  const int64_t _end_ticks;
  double _probability;
  size_t _samples_budget;
  volatile size_t _sample_count;
  Window(double probablity, size_t samples_budget);
 public:
  Window(SamplerWindowParams params, double probability, size_t samples_budget);
  bool should_sample();
  bool should_sample(int64_t timestamp, bool* is_expired);

  size_t sample_count() const {
    const size_t count = Atomic::load(&_sample_count);
    return count < _samples_budget ? count : _samples_budget;
  }

  const SamplerWindowParams params() const {
    return _params;
  }

  /**
   * Ratio between the requested and the measured window duration
   */
   double adjustment_factor() const {
    return (double)(_end_ticks - _start_ticks / now() - _start_ticks);
  }

  double adjustment_factor(jlong window_duration_ms) const;
  bool is_expired() const;

  double probability() const {
    return Atomic::load_acquire(&_probability);
  }

  void set_probability(double probability) {
    Atomic::release_store(&_probability, probability);
  }
};

AdaptiveSampler::Window::Window(SamplerWindowParams params, double probability, size_t samples_budget) :
    _probability(probability),
    _params(params),
    _samples_budget(samples_budget),
    _start_ticks(params.duration == -1 ? 0 : now()),
    _end_ticks(params.duration == -1 ? 0 : _start_ticks + millis_to_countertime(params.duration)),
    _sample_count(0) {}

AdaptiveSampler::Window::Window(double probability, size_t samples_budget) :
  _probability(probability),
  _samples_budget(samples_budget),
  _start_ticks(_params.duration == -1 ? 0 : now()),
  _end_ticks(_params.duration == -1 ? 0 : _start_ticks + millis_to_countertime(_params.duration)),
  _sample_count(0) {}

//   Atomic::inc(&_running_count, memory_order_acq_rel);

bool AdaptiveSampler::Window::should_sample() {
  const double prob = probability();
  const bool sample_all = prob == 1;
  if (sample_all) {
    // if probability is 100% just ignore threshold and always pass
    if (Atomic::add(&_sample_count, (size_t)1, memory_order_acq_rel) <= _samples_budget) {
      return true;
    }
  } else {
    const double n_rand = sample_all ? -1 : get_sampler_support()->next_random_uniform();
    if (n_rand < prob) {
      if (Atomic::add(&_sample_count, (size_t)1, memory_order_acq_rel) <= _samples_budget) {
        return true;
      }
    }
  }
  return false;
}

bool AdaptiveSampler::Window::should_sample(int64_t timestamp, bool* is_expired) {
  assert(is_expired != NULL, "invariant");
  *is_expired = timestamp >= _end_ticks;
  if (*is_expired) {
    return false;
  }
  const double prob = probability();
  const bool sample_all = prob == 1;
  if (sample_all) {
    // if probability is 100% just ignore threshold and always pass
    if (Atomic::add(&_sample_count, (size_t)1, memory_order_acq_rel) <= _samples_budget) {
      return true;
    }
  } else {
    const double n_rand = sample_all ? -1 : get_sampler_support()->next_random_uniform();
    if (n_rand < prob) {
      if (Atomic::add(&_sample_count, (size_t)1, memory_order_acq_rel) <= _samples_budget) {
        return true;
      }
    }
  }
  return false;
}

bool AdaptiveSampler::Window::is_expired() const {
  return now() >= _end_ticks;
}

double AdaptiveSampler::Window::adjustment_factor(jlong window_duration_ms) const {
  return (double)millis_to_countertime(window_duration_ms) / (now() - _start_ticks);
}

static double compute_interval_alpha(size_t interval) {
  return 1 - std::pow(interval, (double)-1 / (double)interval);
}

AdaptiveSampler::AdaptiveSampler(size_t window_lookback_cnt, size_t budget_lookback_cnt) :
    _window_0(NULL),
    _window_1(NULL),
    _active_window(NULL),
    _sampler_support(NULL),
    _window_lookback_alpha(compute_interval_alpha(window_lookback_cnt)),
    _budget_lookback_cnt(budget_lookback_cnt),
    _budget_lookback_alpha(compute_interval_alpha(budget_lookback_cnt)),
    _samples_budget(-1 * (1 + budget_lookback_cnt)),
    _probability(1),
    _avg_samples(std::numeric_limits<double>::quiet_NaN()),
    _avg_count(0),
    _lock(0) {}

bool AdaptiveSampler::initialize() {
  if (_sampler_support == NULL) {
    _sampler_support = new SamplerSupport(true);
    if (_sampler_support == NULL) {
      return false;
    }
  }
  assert(_sampler_support != NULL, "invariant");
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

AdaptiveSampler::Window* AdaptiveSampler::active_window() {
  return Atomic::load_acquire(&_active_window);
}

void AdaptiveSampler::recalculate_averages(SamplerWindowParams new_params) {
  const Window* const window = active_window();
  const SamplerWindowParams previous_params = window->params();
  bool is_dummy = previous_params.duration == -1;
  const double adjustment_factor = is_dummy ? window->adjustment_factor(new_params.duration) : window->adjustment_factor();
  const double samples = window->sample_count() * adjustment_factor;
  // double total_count = _window->total_count() * adjustment_factor;
  double total_count = 0;

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

AdaptiveSampler::Window* AdaptiveSampler::install_new_window(AdaptiveSampler::Window* old) {
  assert(old != NULL, "invariant");
  assert(old == active_window(), "invariant");
  Window* const next = old == _window_0 ? _window_1 : _window_0;
  Atomic::release_store(&_active_window, next);
  return next;
}

void AdaptiveSampler::update_parameters(AdaptiveSampler::Window* window, double probability, size_t samples_budget) {
  assert(window != NULL, "invariant");
  assert(window == active_window(), "invariant");
  window->set_probability(probability);
  window->_samples_budget = samples_budget;
}

// will be called only when the previous window has expired
// exlusively by holder of try_lock
void AdaptiveSampler::rotate_window() {
  Window* const previous = active_window();
  if (!previous->is_expired()) {
    // someone beat us to it
    return;
  }
  Window* const next = install_new_window(previous);
  assert(previous->is_expired(), "invariant");
  SamplerWindowParams params = new_window_params();
  recalculate_averages(params);
  if (_avg_count == 0) {
    _probability = 1;
  } else {
    const double p = (params.sample_count + _samples_budget) / (double)_avg_count;
    _probability = p < 1 ? p : 1;
  }
  update_parameters(next, _probability, _samples_budget);
}

bool AdaptiveSampler::should_sample() {
  bool expired = false;
  const int64_t timestamp = now();
  bool sample = active_window()->should_sample(timestamp, &expired);
  if (expired) {
    {
      JfrTryLock rotation_lock(&_lock);
      if (rotation_lock.acquired()) {
        rotate_window();
      }
    }
    sample = active_window()->should_sample(timestamp, &expired);
  }
  return sample;
}
