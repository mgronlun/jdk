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

class AdaptiveSampler : public JfrCHeapObj {
 private:
  class Window;
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
  void install_new_window(const Window* current_window, SamplerWindowParams next_window_params);
  void rotate_window();

  void recalculate_averages(const Window* current_window, SamplerWindowParams params);
  void recalculate_probability(SamplerWindowParams params);

 protected:
  AdaptiveSampler(size_t window_lookback_cnt, size_t budget_lookback_cnt);
  virtual ~AdaptiveSampler();
  bool initialize();

 public:
  bool should_sample();
  virtual SamplerWindowParams new_window_params() = 0;
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
