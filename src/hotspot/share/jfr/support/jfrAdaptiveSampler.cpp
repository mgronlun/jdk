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
  _start_ticks(0),
  _end_ticks(0),
  _normalization_factor(0),
  _sampling_interval(1),
  _projected_population_size(0),
  _measured_population_size(0) {}

void JfrSamplerWindow::reinitialize(const JfrSamplerParams& params, size_t projected_population_size, size_t sampling_interval) {
  assert(sampling_interval >= 1, "invariant");
  _projected_population_size = projected_population_size;
  Atomic::store(&_measured_population_size, static_cast<size_t>(0));
  if (projected_population_size == 0) {
    // disabled
    _end_ticks = 0;
    return;
  }
  _sampling_interval = sampling_interval;
  _params = params;
  _start_ticks = params.window_duration_ms == 0 ? 0 : now();
  _end_ticks = params.window_duration_ms == 0 ? 0 : _start_ticks + millis_to_countertime(params.window_duration_ms);
}

const JfrSamplerParams& JfrSamplerWindow::params() const {
  return _params;
}

size_t JfrSamplerWindow::population_size_raw() const {
  return Atomic::load(&_measured_population_size);
}

size_t JfrSamplerWindow::population_size() const {
  return floor(_normalization_factor * population_size_raw());
}

size_t JfrSamplerWindow::max_sample_size() const {
  return _projected_population_size / _sampling_interval;
}

size_t JfrSamplerWindow::sample_size(size_t population_size) const {
  if (population_size > _projected_population_size) {
    return max_sample_size();
  }
  return population_size / _sampling_interval;
}

size_t JfrSamplerWindow::sample_size_raw() const {
  return sample_size(population_size_raw());
}

size_t JfrSamplerWindow::sample_size() const {
  return sample_size(population_size());
}

intptr_t JfrSamplerWindow::debt() const {
  return static_cast<intptr_t>(sample_size() - _params.samples_per_window);
}

intptr_t JfrSamplerWindow::cumulative_debt() const {
  return static_cast<intptr_t>((_params.samples_per_window) - max_sample_size()) + debt();
}

double JfrSamplerWindow::duration(int64_t end_ticks) const {
  assert(end_ticks >= _start_ticks, "invariant");
  return static_cast<double>(end_ticks - _start_ticks);
}

JfrAdaptiveSampler::JfrAdaptiveSampler() :
  _window_0(NULL),
  _window_1(NULL),
  _active_window(NULL),
  _support(NULL),
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

inline JfrSamplerWindow* JfrAdaptiveSampler::active_window() const {
  return Atomic::load_acquire(&_active_window);
}

bool JfrAdaptiveSampler::sample(int64_t timestamp) {
  bool expired_window = false;
  const bool result = active_window()->sample(timestamp, &expired_window);
  if (expired_window) {
    {
      JfrTryLock rotate_lock(&_lock);
      if (rotate_lock.acquired()) {
        rotate_window(timestamp);
      }
    }
    return active_window()->sample(now(), &expired_window);
  }
  return result;
}

bool JfrSamplerWindow::sample(int64_t timestamp, bool* expired_window) {
  assert(expired_window != NULL, "invariant");
  assert(*expired_window == false, "invariant");
  if (is_expired(timestamp)) {
    *expired_window = true;
    return false;
  }
  return sample();
}

inline bool JfrSamplerWindow::is_expired(int64_t timestamp) const {
  return timestamp == 0 ? now() >= _end_ticks : timestamp >= _end_ticks;
}

inline bool JfrSamplerWindow::sample() {
  const size_t xchg = Atomic::add(&_measured_population_size, (size_t)1);
  return xchg <= _projected_population_size && xchg % _sampling_interval == 0;
}

// Called exclusively by the holder of the try lock
// when a window is determined to have expired.
void JfrAdaptiveSampler::rotate_window(int64_t timestamp) {
  assert(_lock, "invariant");
  EventSampleWindow event;
  JfrSamplerWindow* const current = active_window();
  if (!current->is_expired(timestamp)) {
    // someone else took care of it
    return;
  }
  current->normalize(timestamp);
  debug(current);
  fill(event, current);
  rotate(current);
  event.commit();
}

inline void JfrSamplerWindow::normalize(int64_t timestamp) {
  _normalization_factor = duration(_end_ticks) / duration(timestamp == 0 ? now() : timestamp);
}

static intptr_t amortize_debt(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  const intptr_t inverse_cumulative_debt = -expired->cumulative_debt(); // negation
  if (params.amortization_windows <= 1) {
    return inverse_cumulative_debt;
  }
  return round(static_cast<double>(inverse_cumulative_debt) / static_cast<double>(params.amortization_windows));
}

static size_t sample_size_for_next_window(const JfrSamplerParams& param, const JfrSamplerWindow* expired) {
  const intptr_t sample_size = static_cast<intptr_t>(param.samples_per_window) + amortize_debt(param, expired);
  return sample_size < 0 ? 0 : sample_size;
}

static size_t sampling_interval_for_next_window(const JfrSamplerWindow* expired, SamplerSupport* support, size_t projected_sample_size) {
  assert(expired != NULL, "invariant");
  assert(support != NULL, "invariant");
  if (projected_sample_size == 0) {
    return 1;
  }
  const size_t previous_population_size = expired->population_size();
  if (previous_population_size < projected_sample_size) {
    return 1;
  }
  const double projected_geometric_mean = static_cast<double>(previous_population_size) / static_cast<double>(projected_sample_size);
  // The sampling interval is a geometric random variable with the specified mean.
  const size_t sampling_interval = support->pick_next_geometric_sample(projected_geometric_mean);
  assert(sampling_interval > 0, "invariant");
  return sampling_interval;
}

void JfrAdaptiveSampler::rotate(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  JfrSamplerParams params = next_window_params(expired); // callback to client for next window parameters
  const size_t projected_sample_size = sample_size_for_next_window(params, expired);
  const size_t sampling_interval = sampling_interval_for_next_window(expired, _support, projected_sample_size);
  assert(sampling_interval >= 1, "invariant");
  const size_t projected_population_size = projected_sample_size * sampling_interval;
  rotate(params, expired, projected_sample_size * sampling_interval, sampling_interval);
}

void JfrAdaptiveSampler::rotate(const JfrSamplerParams& params, const JfrSamplerWindow* expired, size_t projected_population_size, size_t sampling_interval) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  JfrSamplerWindow* const next = next_window(expired);
  assert(next != expired, "invariant");
  next->reinitialize(params, projected_population_size, sampling_interval);
  Atomic::release_store(&_active_window, next);
}

inline JfrSamplerWindow* JfrAdaptiveSampler::next_window(const JfrSamplerWindow* expired) const {
  assert(expired != NULL, "invariant");
  return expired == _window_0 ? _window_1 : _window_0;
}

void JfrAdaptiveSampler::debug(const JfrSamplerWindow* expired) const {
  assert(expired == active_window(), "invariant");
  const JfrSamplerParams& params = expired->params();
  if (params.window_duration_ms == 0) {
    printf("Sampling started...\n");
    return;
  }
  printf("=== sample size: %f, population size: %f, ratio: %f, limit: %f, deviation: %f, cum.deviation: %f\n",
    (double)expired->sample_size(), (double)expired->population_size(),
    expired->population_size() == 0 ? 0 : (double)expired->sample_size() / (double)expired->population_size(),
    (double)expired->max_sample_size(),
    (double)expired->debt(), (double)expired->cumulative_debt());
}

void JfrAdaptiveSampler::fill(EventSampleWindow& event, const JfrSamplerWindow* expired) {
  assert(expired == active_window(), "invariant");
  const JfrSamplerParams& params = expired->params();
  event.set_setPoint(params.samples_per_window);
  event.set_sampleDuration(params.window_duration_ms);
  const size_t sample_size = expired->sample_size();
  event.set_output(sample_size);
  event.set_outputRaw(expired->sample_size_raw());
  const size_t population_size = expired->population_size();
  event.set_input(population_size);
  event.set_ratio(population_size == 0 ? 0 : static_cast<double>(sample_size) / static_cast<double>(population_size));
  event.set_control(0);
  event.set_errorTerm(expired->debt());
  event.set_errors(expired->cumulative_debt());
}
