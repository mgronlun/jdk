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
#include "runtime/atomic.hpp"
#include "runtime/samplerSupport.hpp"
#include "utilities/globalDefinitions.hpp"
#include <cmath>

inline int64_t now() {
  return JfrTicks::now().value();
}

inline double millis_to_countertime(int64_t millis) {
  return static_cast<double>(JfrTimeConverter::nanos_to_countertime(millis * 1000000));
}

JfrSamplerWindow::JfrSamplerWindow() :
  _end_ticks(0),
  _sampling_interval(1),
  _projected_population_size(0),
  _measured_population_size(0) {}

void JfrSamplerWindow::reinitialize(const JfrSamplerParams& params, size_t projected_population_size, size_t sampling_interval) {
  assert(sampling_interval >= 1, "invariant");
  _params = params;
  Atomic::store(&_measured_population_size, static_cast<size_t>(0));
  _projected_population_size = projected_population_size;
  _sampling_interval = sampling_interval;
  const int64_t start_ticks = params.window_duration_ms == 0 ? 0 : now();
  const int64_t end_ticks = start_ticks == 0 ? 0 : start_ticks + millis_to_countertime(params.window_duration_ms);
  Atomic::store(&_end_ticks, end_ticks);
}

size_t JfrSamplerWindow::population_size() const {
  return Atomic::load(&_measured_population_size);
}

size_t JfrSamplerWindow::max_sample_size() const {
  return _projected_population_size / _sampling_interval;
}

size_t JfrSamplerWindow::sample_size() const {
  const size_t size = population_size();
  if (size > _projected_population_size) {
    return max_sample_size();
  }
  return size / _sampling_interval;
}

intptr_t JfrSamplerWindow::debt() const {
  return static_cast<intptr_t>(sample_size() - _params.sample_points_per_window);
}

intptr_t JfrSamplerWindow::accumulated_debt() const {
  return static_cast<intptr_t>((_params.sample_points_per_window) - max_sample_size()) + debt();
}

inline double compute_ewma_alpha_coefficient(size_t lookback_count) {
  return lookback_count <= 1 ? 1 : static_cast<double>(1) / static_cast<double>(lookback_count);
}

JfrAdaptiveSampler::JfrAdaptiveSampler(size_t ewma_lookback_count) :
  _window_0(NULL),
  _window_1(NULL),
  _active_window(NULL),
  _support(NULL),
  _avg_population_size(0),
  _ewma_population_size_alpha(compute_ewma_alpha_coefficient(ewma_lookback_count)),
  _lock(0) {}

JfrAdaptiveSampler::~JfrAdaptiveSampler() {
  delete _window_0;
  delete _window_1;
  delete _support;
}

bool JfrAdaptiveSampler::initialize() {
  assert(_window_0 == NULL, "invariant");
  _window_0 = new JfrSamplerWindow();
  if (_window_0 == NULL) {
    return false;
  }
  assert(_window_1 == NULL, "invariant");
  _window_1 = new JfrSamplerWindow();
  if (_window_1 == NULL) {
    return false;
  }
  _active_window = _window_0;
  assert(_support == NULL, "invariant");
  _support = new SamplerSupport(true);
  return _support != NULL;
}

void JfrAdaptiveSampler::set_window_lookback_count(size_t count) {
  _ewma_population_size_alpha = compute_ewma_alpha_coefficient(count);
}

inline JfrSamplerWindow* JfrAdaptiveSampler::active_window() const {
  return Atomic::load_acquire(&_active_window);
}

bool JfrAdaptiveSampler::sample(int64_t timestamp) {
  bool expired_window;
  const bool result = active_window()->sample(timestamp, &expired_window);
  if (expired_window) {
    JfrTryLock rotate_lock(&_lock);
    if (rotate_lock.acquired()) {
      rotate_window(timestamp);
    }
  }
  return result;
}

bool JfrSamplerWindow::sample(int64_t timestamp, bool* expired_window) {
  assert(expired_window != NULL, "invariant");
  *expired_window = is_expired(timestamp);
  return *expired_window ? false : sample();
}

inline bool JfrSamplerWindow::is_expired(int64_t timestamp) const {
  const int64_t end_ticks = Atomic::load(&_end_ticks);
  return timestamp == 0 ? now() >= end_ticks : timestamp >= end_ticks;
}

inline bool JfrSamplerWindow::sample() {
  const size_t ordinal = Atomic::add(&_measured_population_size, (size_t)1);
  return ordinal <= _projected_population_size && ordinal % _sampling_interval == 0;
}

// Called exclusively by the holder of the try lock
// when a window is determined to have expired.
void JfrAdaptiveSampler::rotate_window(int64_t timestamp) {
  assert(_lock, "invariant");
  EventSamplerWindow event;
  JfrSamplerWindow* const current = active_window();
  if (!current->is_expired(timestamp)) {
    // Someone else took care of it.
    return;
  }
  // debug(current, _avg_population_size);
  fill(event, current);
  rotate(current);
  event.commit();
}

void JfrAdaptiveSampler::debug(const JfrSamplerWindow* expired, double avg_population_size) const {
  assert(expired == active_window(), "invariant");
  const JfrSamplerParams& params = expired->params();
  if (params.window_duration_ms == 0) {
    printf("Sampling started...\n");
    return;
  }
  printf("=== sample size: %f, population size: %f, ratio: %f, limit: %f, deviation: %f, cum.deviation: %f, avg pop size: %f\n",
    (double)expired->sample_size(), (double)expired->population_size(),
    expired->population_size() == 0 ? 0 : (double)expired->sample_size() / (double)expired->population_size(),
    (double)expired->max_sample_size(),
    (double)expired->debt(), (double)expired->accumulated_debt(), avg_population_size);
}

void JfrAdaptiveSampler::fill(EventSamplerWindow& event, const JfrSamplerWindow* expired) {
  assert(expired == active_window(), "invariant");
  const JfrSamplerParams& params = expired->params();
  event.set_setPoint(params.sample_points_per_window);
  event.set_windowDuration(params.window_duration_ms);
  const size_t sample_size = expired->sample_size();
  event.set_sampleSize(sample_size);
  event.set_sampleSizeRaw(expired->sample_size());
  const size_t population_size = expired->population_size();
  event.set_populationSize(population_size);
  event.set_ratio(population_size == 0 ? 0 : static_cast<double>(sample_size) / static_cast<double>(population_size));
  event.set_debt(expired->debt());
  event.set_accumulatedDebt(expired->accumulated_debt());
  event.set_lookbackCount(1 / _ewma_population_size_alpha);
}

static size_t amortization_payment(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  const size_t inverse_accumulated_debt = -expired->accumulated_debt(); // negation
  if (params.amortization_window_count <= 1) {
    return inverse_accumulated_debt;
  }
  return floor(static_cast<double>(inverse_accumulated_debt) / static_cast<double>(params.amortization_window_count));
}

static size_t sample_size_for_next_window(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  return params.sample_points_per_window + amortization_payment(params, expired);
}

static size_t sampling_interval_for_next_window(SamplerSupport* support, size_t projected_sample_size, double avg_population_size) {
  assert(support != NULL, "invariant");
  if (projected_sample_size <= 0) {
    return 1;
  }
  const size_t projected_geometric_mean = avg_population_size / projected_sample_size;
  // The sampling interval is a geometric random variable with the specified mean.
  const size_t sampling_interval = support->pick_next_geometric_sample(projected_geometric_mean);
  // printf("projected geometric mean: %f, sampling interval %llu\n", projected_geometric_mean, sampling_interval);
  assert(sampling_interval > 0, "invariant");
  return sampling_interval;
}

/*
 * Exponentially Weighted Moving Average (EWMA):
 *
 * Y is a datapoint (at time t)
 * S is the current EMWA (at time t-1)
 * alpha represents the degree of weighting decrease, a constant smoothing factor between 0 and 1.
 *
 * A higher alpha discounts older observations faster.
 * Returns the new EWMA for S
*/
inline double exponentially_weighted_moving_average(double Y, double alpha, double S) {
  return alpha * Y + (1 - alpha) * S;
}

static double population_size_for_next_window(const JfrSamplerWindow* expired, double ewma_population_size_alpha, double avg_population_size) {
  assert(expired != NULL, "invariant");
  return exponentially_weighted_moving_average(expired->population_size(), ewma_population_size_alpha, avg_population_size);
}

void JfrAdaptiveSampler::rotate(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  const JfrSamplerParams& params = next_window_params(expired);
  const size_t projected_sample_size = sample_size_for_next_window(params, expired);
  _avg_population_size = population_size_for_next_window(expired, _ewma_population_size_alpha, _avg_population_size);
  const size_t sampling_interval = sampling_interval_for_next_window(_support, projected_sample_size, _avg_population_size);
  assert(sampling_interval >= 1, "invariant");
  rotate(params, expired, projected_sample_size, sampling_interval);
}

void JfrAdaptiveSampler::rotate(const JfrSamplerParams& params, const JfrSamplerWindow* expired, size_t projected_sample_size, size_t sampling_interval) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  JfrSamplerWindow* const next = next_window(expired);
  assert(next != expired, "invariant");
  next->reinitialize(params, projected_sample_size * sampling_interval, sampling_interval);
  Atomic::release_store(&_active_window, next);
}

inline JfrSamplerWindow* JfrAdaptiveSampler::next_window(const JfrSamplerWindow* expired) const {
  assert(expired != NULL, "invariant");
  return expired == _window_0 ? _window_1 : _window_0;
}

