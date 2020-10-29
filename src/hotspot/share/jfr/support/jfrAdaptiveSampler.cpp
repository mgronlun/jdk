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
#include "jfr/jfrEvents.hpp"
#include "jfr/support/jfrAdaptiveSampler.hpp"
#include "jfr/utilities/jfrTime.hpp"
#include "jfr/utilities/jfrTimeConverter.hpp"
#include "jfr/utilities/jfrTryLock.hpp"
#include "logging/log.hpp"
#include "runtime/atomic.hpp"
#include "runtime/interfaceSupport.inline.hpp"
#include "runtime/thread.inline.hpp"
#include "utilities/globalDefinitions.hpp"
#include <cmath>

Window::Window(double probability, size_t samples_budget) :
    _start_ticks(_params.duration == -1 ? 0 : AdaptiveSampler::now()),
    _end_ticks(_params.duration == -1 ? 0 : _start_ticks + AdaptiveSampler::millis_to_countertime(_params.duration)),
    _probability(probability),
    _samples_budget(samples_budget),
    _output_count(0),
    _input_count(0) {}

void Window::reinitialize(double probability, size_t samples_budget, SamplerWindowParams params) {
  _probability = probability;
  Atomic::store(&_output_count, static_cast<size_t>(0));
  Atomic::store(&_input_count, static_cast<size_t>(0));
  _samples_budget = samples_budget;
  _params = params;
  _start_ticks = params.duration == -1 ? 0 : AdaptiveSampler::now(),
  _end_ticks = params.duration == -1 ? 0 : _start_ticks + AdaptiveSampler::millis_to_countertime(params.duration);
}

inline double compute_interval_alpha(size_t interval) {
  return 1 - std::pow(interval, static_cast<double>(-1) / static_cast<double>(interval));
}

AdaptiveSampler::AdaptiveSampler(size_t window_lookback_cnt, size_t budget_lookback_cnt) :
    _window_0(NULL),
    _window_1(NULL),
    _active_window(NULL),
    _window_lookback_alpha(compute_interval_alpha(window_lookback_cnt)),
    _budget_lookback_alpha(compute_interval_alpha(budget_lookback_cnt)),
    _samples_budget(-1 * (1 + budget_lookback_cnt)),
    _probability(1),
    _avg_output(std::numeric_limits<double>::quiet_NaN()),
    _budget_lookback_cnt(budget_lookback_cnt),
    _avg_input(0),
    _lock(0) {}

bool AdaptiveSampler::initialize() {
  assert(_window_0 == NULL, "invariant");
  _window_0 = new Window(_probability, static_cast<size_t>(_samples_budget));
  if (_window_0 == NULL) {
    return false;
  }
  assert(_window_1 == NULL, "invariant");
  _window_1 = new Window(_probability, static_cast<size_t>(_samples_budget));
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

inline int64_t AdaptiveSampler::now() {
  return JfrTicks::now().value();
}

inline int64_t AdaptiveSampler::millis_to_countertime(int64_t millis) {
  return JfrTimeConverter::nanos_to_countertime(millis * 1000000);
}

inline Window* AdaptiveSampler::active_window() const {
  return Atomic::load_acquire(&_active_window);
}

bool AdaptiveSampler::should_sample() {
  bool expired_window = false;
  int64_t timestamp = AdaptiveSampler::now();
  bool sample = active_window()->should_sample(timestamp, &expired_window);
  if (expired_window) {
    {
      JfrTryLock rotation_lock(&_lock);
      if (rotation_lock.acquired()) {
        rotate_window();
      }
    }
    timestamp = AdaptiveSampler::now();
    sample = active_window()->should_sample(timestamp, &expired_window);
  }
  return sample;
}

bool Window::should_sample(int64_t timestamp, bool* expired_window) {
  assert(expired_window != NULL, "invariant");
  *expired_window = is_expired(timestamp);
  return *expired_window ? false : enter();
}

bool Window::enter() {
  // Increment number of attempts
  Atomic::add(&_input_count, static_cast<size_t>(1), memory_order_acq_rel);
  return random_trial();
}

/*
 * A window is a fixed time duration together with a target sampling rate / probability
 * value set by the adaptive sampler from what it has learned from the past.
 * The actual sampling rate for the window is modelled using a random process, where an
 * independent Bernoulli trial is performed for each attempt.
 * A trial picks a single uniformly distributed random number in the interval (0,1)
 * for the continuous random variable X. The probability value set for the window
 * is the p value, i.e. the probability of success, for the trial.
 *
 *  If X <= p, the outcome is successful and a sample is taken (if budget permits).
 */

inline double next_random_uniform() {
  return Thread::current()->jfr_thread_local()->sampler_support()->next_random_uniform();
}

bool Window::random_trial() {
  const double p = probability();
  assert(p >= 0, "invariant");
  assert(p <= 1, "invariant");
  if (p < 1) {
    const double X = next_random_uniform();
    if (X > p) {
      // Event did not occur. Not selected for sampling.
      return false;
    }
  }
  // Event occurred. If still within budget, the attempt is selected for sampling.
  return Atomic::add(&_output_count, static_cast<size_t>(1), memory_order_acq_rel) <= _samples_budget;
}

// Called when a window was determined to have expired
// and only exlusively by the holder of the try_lock
void AdaptiveSampler::rotate_window() {
  assert(_lock, "invariant");
  Window* const current_window = active_window();
  if (!current_window->is_expired()) {
    // someone took care of it
    return;
  }
  const SamplerWindowParams new_params = new_window_params();
  EventRetiredSampleWindow event;
  recalculate_averages(current_window, new_params, event);
  recalculate_probability(new_params, event);
  install_next_window(current_window, new_params);
  event.commit();
}

/*
* EWMA:
* 
* Y is a datapoint (at time t)
* S is the current EMWA value (at time t-1)
* alpha represents the degree of weighting decrease, a constant smoothing factor between 0 and 1.
* A higher alpha discounts older observations faster.
*
* returns a new EWMA value for S
*/

inline double exponentially_weighted_moving_average(double Y, double alpha, double S) {
  return std::isnan(S) ? Y : alpha * Y + (1 - alpha) * S;
}

inline double derivative_exponentially_weighted_moving_average(double Y, double alpha, double S) {
  return std::isnan(S) ? Y : S + alpha * (Y - S);
}

void AdaptiveSampler::recalculate_averages(const Window* current_window, SamplerWindowParams new_params, EventRetiredSampleWindow& event) {
  assert(current_window != NULL, "invariant");
  assert(current_window == active_window(), "invariant");
  const SamplerWindowParams current_params = current_window->params();
  const bool is_dummy = current_params.duration == -1;
  const double adjustment_factor = is_dummy ? current_window->adjustment_factor(new_params.duration) : current_window->adjustment_factor();
  const double samples_raw = current_window->output_count();
  const double samples = samples_raw * adjustment_factor;
  const double total_count_raw = current_window->input_count();
  const double total_count = total_count_raw * adjustment_factor;

  const double requested_ratio = (double)current_params.sample_count / (double)current_params.duration;

  event.set_requested_ratio(requested_ratio);
  event.set_samples_raw(samples_raw);
  event.set_attempts_raw(total_count_raw);
  event.set_raw_ratio(samples_raw / total_count_raw);
  event.set_samples(samples);
  event.set_attempts(total_count);
  event.set_adjustment_factor(adjustment_factor);
  event.set_adjusted_ratio(samples / total_count);
  event.set_prev_probability(current_window->probability());
  event.set_prev_avgSamples(_avg_output);
  event.set_prev_avgAttempts(_avg_input);
  event.set_prev_budget(current_window->_samples_budget);

  if (!is_dummy) {
    _avg_output = exponentially_weighted_moving_average(samples, _budget_lookback_alpha, _avg_output);
  }
  _samples_budget = fmax<double>(new_params.sample_count - _avg_output, 0) * _budget_lookback_cnt;
  printf("=== avg samples: %f, samples: %f, prob: %f\n", _avg_output, samples, current_window->probability());
  _avg_input = _avg_input == 0 ? total_count : exponentially_weighted_moving_average(total_count, _window_lookback_alpha, static_cast<double>(_avg_input));

  /*

  if (!is_dummy) {
    _avg_output = std::isnan(_avg_output) ? samples : _avg_output + _budget_lookback_alpha * (samples - _avg_output);
  }
  _samples_budget = fmax<double>(new_params.sample_count - _avg_output, 0) * _budget_lookback_cnt;

  printf("=== samples raw: %f, attempts_raw: %f, ratio: %f, prob: %f\n", samples_raw, total_count_raw, samples_raw / total_count_raw, current_window->probability());
  printf("=== avg samples: %f, samples: %f, adjustment: %f\n", _avg_output, samples, adjustment_factor);

  if (_avg_input == 0) {
    _avg_input = total_count;
  } else {
    // need to convert int '*_count' variables to double to prevent bit overflow
    _avg_input = _avg_input + _window_lookback_alpha * ((double)total_count - (double)_avg_input);
  }

  */

  event.set_new_avgSamples(_avg_output);
  event.set_new_avgAttempts(_avg_input);

  /*
  const double samples = current_window->output_count() * adjustment_factor;
  const double attempts = current_window->input_count() * adjustment_factor;
  if (!is_dummy) {
    _avg_output = derivative_exponentially_weighted_moving_average(samples, _budget_lookback_alpha, _avg_output);
  }
  _samples_budget = fmax<double>(new_params.sample_count - _avg_output, 0) * _budget_lookback_cnt;
  printf("=== avg samples: %f, samples: %f, adjustment: %f\n", _avg_output, samples, adjustment_factor);
  _avg_input = _avg_input == 0 ? attempts : derivative_exponentially_weighted_moving_average(attempts, _window_lookback_alpha, static_cast<double>(_avg_input));
  */
}

/*
double samples = current_window->output_count() * adjustment_factor;
double total_count = current_window->input_count() * adjustment_factor;

if (!is_dummy) {
  _avg_output = std::isnan(_avg_output) ? samples : _avg_output + _budget_lookback_alpha * (samples - _avg_output);
}
_samples_budget = fmax<double>(new_params.sample_count - _avg_output, 0) * _budget_lookback_cnt;

printf("=== avg samples: %f, samples: %f, adjustment: %f\n", _avg_output, samples, adjustment_factor);

if (_avg_input == 0) {
  _avg_input = total_count;
}
else {
  // need to convert int '*_count' variables to double to prevent bit overflow
  _avg_input = _avg_input + _window_lookback_alpha * ((double)total_count - (double)_avg_input);
}
*/

void AdaptiveSampler::recalculate_probability(SamplerWindowParams params, EventRetiredSampleWindow& event) {
  if (_avg_input == 0) {
    _probability = 1;
  } else {
    const double p = (params.sample_count + _samples_budget) / static_cast<double>(_avg_input);
    _probability = p < 1 ? p : 1;
  }
  event.set_new_probability(_probability);
  event.set_new_budget(_samples_budget);
}

inline Window* AdaptiveSampler::next_window(const Window* current_window) const {
  assert(current_window != NULL, "invariant");
  return current_window == _window_0 ? _window_1 : _window_0;
}

void AdaptiveSampler::install_next_window(const Window* current_window, SamplerWindowParams new_params) {
  assert(current_window != NULL, "invariant");
  assert(current_window == active_window(), "invariant");
  Window* const next = next_window(current_window);
  assert(next != current_window, "invariant");
  next->reinitialize(_probability, _samples_budget, new_params);
  Atomic::release_store(&_active_window, next);
}

FixedRateSampler::FixedRateSampler(int64_t window_duration, int64_t samples_per_window, size_t window_lookback_cnt, size_t budget_lookback_cnt) :
    AdaptiveSampler(window_lookback_cnt, budget_lookback_cnt) {
  _params.sample_count = samples_per_window;
  _params.duration = window_duration;
}

