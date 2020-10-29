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

#ifndef SHARE_JFR_SUPPORT_JFRADAPTIVESAMPLER_HPP
#define SHARE_JFR_SUPPORT_JFRADAPTIVESAMPLER_HPP

#include "jfr/utilities/jfrAllocation.hpp"
#include "runtime/atomic.hpp"

/**
 * The adaptive sampler is guaranteeing a maximum number of samples picked per a certain time interval.
 * The maximum number is a soft limit which can be, in extreme situations be crossed but the overshoot
 * is usually staying within 15-20% of the requested limit (the actual number is affected by the 'sampling budget' size).
 * As for the implementation - the sampler is using fixed size time windows and adjusts the sampling rate/probability
 * for the next window based on what it learned in the past. While this strategy fares pretty well for a fairly stable system
 * it can fail for bursty ones when an extremely bursty window can influence the moving average in a way that several subsequent
 * windows will end up undersampled. As a measure of compensation the adaptive sampler employs the concept of 'sampling budget'
 * The 'sampling budget' is working as a 'spike damper', smoothing out the extremes in a way that the overall target rate
 * is obeyed without highly over- or under-sampled winows.
 */

struct SamplerWindowParams {
  int64_t duration;
  int64_t sample_count;
};

class EventRetiredSampleWindow;
class Window;

class AdaptiveSampler : public JfrCHeapObj {
 friend class Window;
 private:
  Window* _window_0;
  Window* _window_1;
  Window* _active_window;
  const double _window_lookback_alpha;
  const double _budget_lookback_alpha;
  double _samples_budget;
  double _probability;
  double _avg_output;
  const size_t _budget_lookback_cnt;
  size_t _avg_input;
  volatile int _lock;

  Window* active_window() const;
  Window* next_window(const Window* current_window) const;
  void install_next_window(const Window* current_window, SamplerWindowParams next_window_params);
  void rotate_window();

  void recalculate_averages(const Window* current_window, SamplerWindowParams params, EventRetiredSampleWindow& event);
  void recalculate_probability(SamplerWindowParams params, EventRetiredSampleWindow& event);

  static int64_t now();
  static int64_t millis_to_countertime(int64_t millis);

 protected:
  AdaptiveSampler(size_t window_lookback_cnt, size_t budget_lookback_cnt);
  virtual ~AdaptiveSampler();

 public:
  bool initialize();
  bool should_sample();
  virtual SamplerWindowParams new_window_params() = 0;
};

class Window : public JfrCHeapObj {
  friend class AdaptiveSampler;
 private:
  SamplerWindowParams _params = { -1, -1 };
  int64_t _start_ticks;
  int64_t _end_ticks;
  double _probability;
  size_t _samples_budget;
  volatile size_t _output_count;
  volatile size_t _input_count;

  Window(double probablity, size_t samples_budget);
  void reinitialize(double probability, size_t samples_budget, SamplerWindowParams params);

  bool should_sample(int64_t timestamp, bool* is_expired);
  bool enter();
  bool random_trial();

  inline
  size_t output_count() const {
    const size_t count = Atomic::load(&_output_count);
    return count < _samples_budget ? count : _samples_budget;
  }

  inline
  size_t input_count() const {
    return Atomic::load(&_input_count);
  }

  inline
  const SamplerWindowParams params() const {
    return _params;
  }

  double probability() const {
    return Atomic::load_acquire(&_probability);
  }

  /**
   * Ratio between the requested and the measured window duration
   */
  inline
  double adjustment_factor() const {
    return static_cast<double>(_end_ticks - _start_ticks) / static_cast<double>(AdaptiveSampler::now() - _start_ticks);
  }

  inline
  double adjustment_factor(int64_t window_duration_ms) const {
    return static_cast<double>(AdaptiveSampler::millis_to_countertime(window_duration_ms)) / static_cast<double>(AdaptiveSampler::now() - _start_ticks);
  }

  inline
  bool is_expired() const {
    return is_expired(AdaptiveSampler::now());
  }

  inline
  bool is_expired(int64_t timestamp) const {
    return timestamp >= _end_ticks;
  }
};

class FixedRateSampler : public AdaptiveSampler {
 private:
  SamplerWindowParams _params;
 public:
  FixedRateSampler(int64_t window_duration, int64_t samples_per_window, size_t window_lookback_cnt, size_t budget_lookback_cnt);
  SamplerWindowParams new_window_params() {
    return _params;
  }
};

#endif // SHARE_JFR_SUPPORT_JFRADAPTIVESAMPLER_HPP
