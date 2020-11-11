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

/*
 * The terminology is mostly from the domain of statistics:
 *
 * Population - a set of elements of interest.
 * Sample - a subset of elements from a population selected by a defined procedure.
 * Sample point - an element of a sample (sub)set.
 * Sampling interval - the distance or time between which measurements are taken, also referred to as 'nth selection'
 * Debt - an error term, signifying the deviation from a configured set point.
 * Amortization - a projection or strategy to recover accumulated debt.
 * Window - as in time window or time frame. The sampler sees the evolution of the system in time slices, i.e. in windows.
 * Rotate - the process of retiring an expired window and installing a new window. Similar somewhat to context switching.
 *
 * The adaptive sampler will, on average, guarantee a maximum number of sample points selected from a populuation
 * over a certain time interval. It is using fixed size time windows and adjusts the sampling interval for the next
 * window based on what it learned in the past. Each window has a set point, which is the target number of sample points
 * to select on average. The sampler keeps a cumulative error term, also called 'accumulated debt', which is a measure
 * for how much the sampler is deviating from the set point over time. The maximum number of sample points selected
 * during an individual window is the set point + the accumulated debt.
 * Hence, the 'accumulated debt' is also working as a 'spike damper', smoothing out the extremes in a way that the overall
 * target rate is obeyed without highly over- or under-sampled windows.
 *
 * Sample point selection is defined by a sampling interval, which gives instructions for selecting the 'nth' element
 * in a population. The sampling interval is a geometric random variable, determined by a stochastic process and
 * recalculated for each window in an effort keep the sample set representative.
 *
 * Each window is configured individually, by an instance of the SamplerParams struct. On window expiration,
 * but before switching in the next window, the sampler calls a subclass with the just expired window as an argument.
.* A subclass can inspect the window to study the history of the system and also get
 * an overview of how the sampler is performing to help draw inferences. Based on what it learned, it can choose to
 * let the sampler re-apply an updated set of parameters to the next, upcoming, window. This is a basic feedback control loop
 * that can be developed further, perhaps evolving more elaborate sampling schemes in the future.
 *
 * Using the AdaptiveSampler mechanism, a user can now specify at a high level, for example that he/she would like a
 * maximum rate of n sample points per second. Note that the sampler only guarantees a maxmimum rate of n on average.
 * Naturally, lower rates will be reported if the system does not produce a population to sustain the requested rate,
 * but n per second is respected as a maximum limit hence it will never report an average rate higher than n per second.
 *
 * One good use of the sampler is to employ it as a throttler, or regulator, to help shape large data sets into smaller,
 * more managable subsets while still keeping the data somewhat representative.
 *
 */

struct JfrSamplerParams {
  size_t sample_points_per_window;  // the number of sample points to select per window on average, the set point for the sampler
  size_t window_duration_ms;        // period of time the sampler selects n sample points from a population
  size_t amortization_window_count; // future projection / strategy to recover accumulated debt (specified in number of windows)
};

class JfrSamplerWindow : public JfrCHeapObj {
  friend class JfrAdaptiveSampler;
 private:
  JfrSamplerParams _params = { 0, 0, 0 };
  volatile int64_t _end_ticks;
  size_t _sampling_interval;
  size_t _projected_population_size;
  volatile size_t _measured_population_size;

  JfrSamplerWindow();
  void reinitialize(const JfrSamplerParams& params, size_t projected_population_size, size_t sampling_interval);
  bool is_expired(int64_t timestamp) const;
  size_t max_sample_size() const;
  bool sample(int64_t timestamp, bool* is_expired);
  bool sample();

 public:
  size_t population_size() const;
  size_t sample_size() const;
  intptr_t debt() const;
  intptr_t accumulated_debt() const;
  const JfrSamplerParams& params() const {
    return _params;
  }
};

class EventSamplerWindow;
class SamplerSupport;

class JfrAdaptiveSampler : public JfrCHeapObj {
 private:
  JfrSamplerWindow* _window_0;
  JfrSamplerWindow* _window_1;
  JfrSamplerWindow* _active_window;
  SamplerSupport* _support;
  double _avg_population_size;
  double _ewma_population_size_alpha;
  volatile int _lock;
  bool _reset;
 protected:
  JfrAdaptiveSampler(size_t window_lookback_count);
  virtual ~JfrAdaptiveSampler();
  virtual bool initialize();
  void reset();

  JfrSamplerWindow* active_window() const;
  JfrSamplerWindow* next_window(const JfrSamplerWindow* expired) const;
  void set_window_lookback_count(size_t count);
  void rotate_window(int64_t timestamp);
  void rotate(const JfrSamplerWindow* expired);
  void rotate(const JfrSamplerParams& params, const JfrSamplerWindow* expired, size_t projected_sample_size, size_t sampling_interval);
  virtual const JfrSamplerParams& next_window_params(const JfrSamplerWindow* expired) = 0;

  void debug(const JfrSamplerWindow* expired, double avg_population_size) const;
  void fill(EventSamplerWindow& event, const JfrSamplerWindow* expired);

 public:
  bool sample(int64_t timestamp = 0);
};

class JfrFixedRateSampler : public JfrAdaptiveSampler {
 private:
  JfrSamplerParams _params;
 public:
  JfrFixedRateSampler(size_t samples_per_window, size_t window_duration_ms);
  const JfrSamplerParams& next_window_params(const JfrSamplerWindow* expired) {
    return _params;
  }
};

#endif // SHARE_JFR_SUPPORT_JFRADAPTIVESAMPLER_HPP
