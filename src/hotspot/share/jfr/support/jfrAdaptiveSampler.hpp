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
 * is obeyed without highly over- or under-sampled windows.
 */

struct JfrSamplerParams {
  size_t samples_per_window;   // the average number of sample points per window
  size_t window_duration_ms;   // period of time the sampler selects "samples_per_window" number of samples
  size_t amortization_windows; // future projection / strategy to recover the cumulative deviation (specified in number of windows)
};

class JfrSamplerWindow : public JfrCHeapObj {
  friend class JfrAdaptiveSampler;
 private:
  JfrSamplerParams _params = { 0, 0, 0 };
  int64_t _start_ticks;
  int64_t _end_ticks;
  size_t _sampling_interval;
  size_t _projected_population_size;
  volatile size_t _measured_population_size;
  double _normalization_factor;

  size_t population_size_raw() const;
  size_t sample_size_raw() const;
  size_t sample_size(size_t population_size) const;
  size_t max_sample_size() const;

  JfrSamplerWindow();
  void reinitialize(const JfrSamplerParams& params, size_t projected_population_size, size_t sampling_interval);
  double duration(int64_t end_ticks) const;
  bool is_expired(int64_t timestamp) const;
  void normalize(int64_t timestamp);

  bool sample(int64_t timestamp, bool* is_expired);
  bool sample();

 public:
  const JfrSamplerParams& params() const;
  size_t population_size() const;
  size_t sample_size() const;
  intptr_t debt() const;
  intptr_t cumulative_debt() const;
};

class EventSampleWindow;
class SamplerSupport;

class JfrAdaptiveSampler : public JfrCHeapObj {
 protected:
  JfrSamplerWindow* _window_0;
  JfrSamplerWindow* _window_1;
  JfrSamplerWindow* _active_window;
  SamplerSupport* _support;
  volatile int _lock;

  JfrAdaptiveSampler();
  virtual ~JfrAdaptiveSampler();
  virtual bool initialize();

  void debug(const JfrSamplerWindow* expired) const;
  void fill(EventSampleWindow& event, const JfrSamplerWindow* expired);

  JfrSamplerWindow* active_window() const;
  JfrSamplerWindow* next_window(const JfrSamplerWindow* expired) const;
  void rotate_window(int64_t timestamp);
  void rotate(const JfrSamplerWindow* expired);
  void rotate(const JfrSamplerParams& params, const JfrSamplerWindow* expired, size_t sampling_interval, size_t projected_population_size);

  virtual JfrSamplerParams next_window_params(const JfrSamplerWindow* expired) = 0;

 public:
  virtual bool sample(int64_t timestamp = 0);
};

class JfrFixedRateSampler : public JfrAdaptiveSampler {
 private:
  JfrSamplerParams _params;
 public:
  JfrFixedRateSampler(size_t samples_per_window, size_t window_duration_ms);
  JfrSamplerParams next_window_params(const JfrSamplerWindow* expired) {
    return _params;
  }
};

#endif // SHARE_JFR_SUPPORT_JFRADAPTIVESAMPLER_HPP
