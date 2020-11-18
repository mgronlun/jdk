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
#include "runtime/thread.inline.hpp"
#include "utilities/globalDefinitions.hpp"
#include <cmath>

constexpr const double MAX_PROBABILITY = 1.0;
constexpr const size_t MAX_SIZE_T = static_cast<size_t>(-1);

inline double next_random_uniform() {
  return Thread::current()->jfr_thread_local()->sampler_support()->next_random_uniform();
}

class JfrProbabilityWindow : public JfrSamplerWindow {
  friend class JfrAdaptiveSampler;
 private:
  double _probability;
  size_t _max_sample_size;

  JfrProbabilityWindow() : JfrSamplerWindow(), _probability(MAX_PROBABILITY), _max_sample_size(0) {}

  void reset() override {
    _probability = MAX_PROBABILITY;
    JfrSamplerWindow::reset();
  }

  bool is_derived() const override {
    return true;
  }

  double probability() const {
    return Atomic::load(&_probability);
  }

  void set_probability(double p) {
    assert(p >= 0.0, "invariant");
    assert(p <= MAX_PROBABILITY, "invariant");
    Atomic::store(&_probability, p);
  }

  /*
   * A Bernoulli trial assigns a single uniformly distributed random number
   * in the interval [0,1] to the continuous random variable X. The probability
   * value assigned to the window is the p value of the trial,
   * i.e. the probability of success.
   *
   *  If X < p, the outcome is successful.
   */
  bool sample() const override {
    const double p = probability();
    if (p == MAX_PROBABILITY) return true;
    if (measured_sample_size() >= _max_sample_size) {
      return false;
    }
    if (next_random_uniform() >= p) {
      // Event did not occur.
      return false;
    }
    return Atomic::add(&_measured_sample_size, static_cast<size_t>(1)) <= _max_sample_size;
  }

 public:
  size_t sample_size() const override {
    const size_t measured_size = measured_sample_size();
    return measured_size < _max_sample_size ? measured_size : _max_sample_size;
  }
};

inline int64_t now() {
  return JfrTicks::now().value();
}

inline int64_t millis_to_countertime(int64_t millis) {
  return JfrTimeConverter::nanos_to_countertime(millis * NANOSECS_PER_MILLISEC);
}

JfrSamplerWindow::JfrSamplerWindow() :
  _params(),
  _end_ticks(0),
  _nth_mod_value(0),
  _sampling_interval(1),
  _projected_sample_size(0),
  _projected_population_size(0),
  _measured_sample_size(0),
  _measured_population_size(0) {
  _params.sample_points_per_window = 0;
  _params.probability = JfrSamplerParams::unused;
  _params.nth_selection = JfrSamplerParams::unused;
}

void JfrSamplerWindow::initialize(const JfrSamplerParams& params, size_t initial_value) {
  _params = params;
  if (params.window_duration_ms == 0) {
    Atomic::store(&_end_ticks, static_cast<int64_t>(0));
    return;
  }
  Atomic::store(&_measured_sample_size, static_cast<size_t>(0));
  Atomic::store(&_measured_population_size, initial_value);
  const int64_t end_ticks = now() + millis_to_countertime(params.window_duration_ms);
  Atomic::store(&_end_ticks, end_ticks);
}

void JfrSamplerWindow::reset() {
  _nth_mod_value = 0;
  _sampling_interval = 1;
  _projected_sample_size = MAX_SIZE_T;
  _projected_population_size = MAX_SIZE_T;
}

inline bool JfrSamplerWindow::is_derived() const {
  return false;
}

size_t JfrSamplerWindow::population_size() const {
  return Atomic::load(&_measured_population_size);
}

size_t JfrSamplerWindow::measured_sample_size() const {
  return Atomic::load(&_measured_sample_size);
}

size_t JfrSamplerWindow::max_sample_size() const {
  return _projected_sample_size;
}

size_t JfrSamplerWindow::sample_size() const {
  // For nth selection and discrete probabilities, the sample size
  // is derived from the measured population size.
  const size_t size = population_size();
  if (size > _projected_population_size) {
    return max_sample_size();
  }
  if (size == 0 || size < _nth_mod_value) {
    return 0;
  }
  return _nth_mod_value == 0 ? size / _sampling_interval : ((size - _nth_mod_value) / _sampling_interval) + 1;
}

intptr_t JfrSamplerWindow::debt() const {
  return static_cast<intptr_t>(sample_size() - _params.sample_points_per_window);
}

intptr_t JfrSamplerWindow::accumulated_debt() const {
  return static_cast<intptr_t>(_params.sample_points_per_window - max_sample_size()) + debt();
}

JfrAdaptiveSampler::JfrAdaptiveSampler() :
  _window_0(NULL),
  _window_1(NULL),
  _window_2(NULL),
  _window_3(NULL),
  _active_window(NULL),
  _avg_population_size(0),
  _ewma_population_size_alpha(0),
  _acc_debt_carry_limit(0),
  _acc_debt_carry_count(0),
  _initial_value_for_next_window(0),
  _lock(0),
  _probability_mode(false) {}

JfrAdaptiveSampler::~JfrAdaptiveSampler() {
  delete _window_0;
  delete _window_1;
  delete _window_2;
  delete _window_3;
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
  // Specialized windows for handling continuous probabilities
  assert(_window_2 == NULL, "invariant");
  _window_2 = new JfrProbabilityWindow();
  if (_window_2 == NULL) {
    return false;
  }
  assert(_window_3 == NULL, "invariant");
  _window_3 = new JfrProbabilityWindow();
  if (_window_3 == NULL) {
    return false;
  }
  _active_window = _window_0;
  return true;
}

inline const JfrSamplerWindow* JfrAdaptiveSampler::active_window() const {
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

bool JfrSamplerWindow::sample(int64_t timestamp, bool* expired_window) const {
  assert(expired_window != NULL, "invariant");
  *expired_window = is_expired(timestamp);
  return *expired_window ? false : sample();
}

inline bool JfrSamplerWindow::is_expired(int64_t timestamp) const {
  const int64_t end_ticks = Atomic::load(&_end_ticks);
  return timestamp == 0 ? now() >= end_ticks : timestamp >= end_ticks;
}

inline bool JfrSamplerWindow::sample() const {
  const size_t ordinal = Atomic::add(&_measured_population_size, static_cast<size_t>(1));
  return ordinal <= _projected_population_size && ordinal % _sampling_interval == _nth_mod_value;
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

static double next_window_population_size(const JfrSamplerWindow* expired, double alpha, double avg_population_size) {
  assert(expired != NULL, "invariant");
  return exponentially_weighted_moving_average(expired->population_size(), alpha, avg_population_size);
}

/*
 * Based on previous information, the sampler creates a future 'projection',
 * a speculation of what the situation will be for the next window.
 * Given this projection, parameters are set accordingly to collect
 * a sample set as close as possible to the target (set point), which
 * is a function of the number of sample_points_per_window + accumulated debt.
 */
void JfrAdaptiveSampler::set_rate(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(!_probability_mode, "invariant");
  const size_t projected_sample_size = next_window_sample_size(params, expired);
  if (projected_sample_size == 0) {
    return;
  }
  JfrSamplerWindow* const next = next_window(expired);
  assert(next != expired, "invariant");
  // _avg_population_size = next_window_population_size(expired, _ewma_population_size_alpha, _avg_population_size);
  if (expired->population_size() > projected_sample_size) {
    next->_sampling_interval = expired->population_size() / projected_sample_size;
    // next->_nth_mod_value = randomized_nth_selection(next->_sampling_interval);
  } else {
    next->_sampling_interval = 1;
    //next->_nth_mod_value = 0;
  }
  next->_nth_mod_value = 0;
  // set_projected_population_size(projected_sample_size, next);
  next->_projected_sample_size = projected_sample_size;
  next->_projected_population_size = projected_sample_size * next->_sampling_interval;
}

size_t JfrAdaptiveSampler::next_window_sample_size(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::JfrSamplerParams::unused, "invariant");
  if (params.sample_points_per_window == 0) {
    next_window(expired)->_projected_population_size = 0;
    return 0;
  }
  return params.sample_points_per_window + next_window_amortization(expired);
}

void JfrAdaptiveSampler::set_projected_population_size(size_t projected_sample_size, JfrSamplerWindow* next) {
  assert(projected_sample_size > 0, "invariant");
  assert(next != NULL, "invariant");
  assert(next != active_window(), "invariant");
  assert(next->_sampling_interval > 0, "invariant");
  assert(next->_nth_mod_value < next->_sampling_interval, "invariant");
  next->_projected_sample_size = projected_sample_size;
  if (next->_nth_mod_value == 0) {
    next->_projected_population_size = projected_sample_size * next->_sampling_interval;
    return;
  }
  next->_projected_population_size = ((projected_sample_size - 1) * next->_sampling_interval) + next->_nth_mod_value;
}

/*
 * When the sampler is configured to maintain a rate, is employs the concepts
 * of 'debt' and 'accumulated debt'. 'Accumulated debt' can be thought of as
 * a cumulative error term, and is indicative for how much the sampler is
 * deviating from a set point (the rate). Debt accumulates naturally over time,
 * as a function of undersampled windows, mainly because of system fluctuations,
 * i.e. too small populations.
 *
 * A specified rate is implicitly a _maximal_ rate, so the sampler must ensure
 * to respect this 'limit'. Rates are normalized to per-second ratios, so the
 * 'limit' to respect is on a per second basis.During this second, the sampler
 * has freedom to dynamically re-adjust, and it does so by 'amortizing'
 * accumulated debt over a certain number of windows that fall within the second.
 *
 * Intuitively, accumulated debt 'carry over' from the predecessor to the successor
 * window iff within the allowable time frame (determined by _acc_debt_carry_limit).
 * The successor window will sample more points to make amends, or 'amortize' debt
 * accumulated by its predecessor(s).
 */
size_t JfrAdaptiveSampler::next_window_amortization(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  const size_t accumulated_debt = -expired->accumulated_debt(); // negation
  return accumulated_debt;
  /*
  if (_acc_debt_carry_count == _acc_debt_carry_limit) {
    _acc_debt_carry_count = 1;
    return 0;
  }
  ++_acc_debt_carry_count;
  const size_t accumulated_debt = -expired->accumulated_debt(); // negation
  return accumulated_debt;
  /*
  if (_acc_debt_carry_count == _acc_debt_carry_limit) {
    return accumulated_debt;
  }
  return accumulated_debt / (_acc_debt_carry_limit - _acc_debt_carry_count);
  */
}

inline bool is_discrete_probability(double p, size_t* nth_selection) {
  assert(nth_selection != NULL, "invariant");
  assert(p >= 0.0, "invariant");
  assert(p <= MAX_PROBABILITY, "invariant");
  if (p == 0.0) return false;
  const double p_interval = static_cast<double>(1) / p;
  *nth_selection = p_interval;
  return p_interval == static_cast<double>(*nth_selection);
}

// Randomly select the nth element in the interval, map it to its zero based index.
inline size_t randomized_nth_selection_mod_value(size_t interval) {
  assert(interval > 0, "invariant");
  if (interval == 1) {
    return 0;
  }
  const double factor = next_random_uniform();
  assert(factor <= 1.0, "invariant");
  return factor * (interval - 1);
}

void JfrAdaptiveSampler::set_rate_and_probability(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  const size_t projected_sample_size = next_window_sample_size(params, expired);
  if (projected_sample_size == 0) {
    return;
  }
  size_t probability_nth_selection;
  if (is_discrete_probability(params.probability, &probability_nth_selection)) {
    assert(probability_nth_selection > 0, "invariant");
    // For discrete probabilities, we can use randomized nth selection.
    JfrSamplerWindow* const next = next_window(expired);
    next->_sampling_interval = probability_nth_selection;
    next->_nth_mod_value = randomized_nth_selection_mod_value(probability_nth_selection);
    set_projected_population_size(projected_sample_size, next);
    return;
  }
  set_probability(params, projected_sample_size, expired);
}

void JfrAdaptiveSampler::set_rate_and_nth_selection(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  const size_t projected_sample_size = next_window_sample_size(params, expired);
  if (projected_sample_size == 0) {
    return;
  }
  JfrSamplerWindow* const next = next_window(expired);
  next->_sampling_interval = params.nth_selection;
  set_projected_population_size(projected_sample_size, next);
}

void JfrAdaptiveSampler::set_nth_selection(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  next_window(expired)->_sampling_interval = params.nth_selection;
  // To keep a consistent modulo over the window boundary.
  _initial_value_for_next_window = expired->_sampling_interval == 0 ? 0 :
                                     expired->population_size() % expired->_sampling_interval;
}

// Continuous probabilities employ a Bernoulli trial.
void JfrAdaptiveSampler::set_probability(const JfrSamplerParams& params, size_t projected_sample_size, const JfrSamplerWindow* expired) {
  _probability_mode = true;
  JfrProbabilityWindow* const next_prob_window = static_cast<JfrProbabilityWindow*>(next_window(expired));
  next_prob_window->_max_sample_size = projected_sample_size;
  next_prob_window->set_probability(params.probability);
}

void JfrAdaptiveSampler::set_probability(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  size_t probability_nth_selection;
  if (is_discrete_probability(params.probability, &probability_nth_selection)) {
    assert(probability_nth_selection > 0, "invariant");
    JfrSamplerWindow* const next = next_window(expired);
    // For discrete probabilities, we can use randomized nth selection for better performance.
    next->_sampling_interval = probability_nth_selection;
    next->_nth_mod_value = randomized_nth_selection_mod_value(probability_nth_selection);
    assert(next->_nth_mod_value < next->_sampling_interval, "invariant");
    return;
  }
  set_probability(params, params.probability == 0 ? 0 : MAX_SIZE_T, expired);
}

void JfrAdaptiveSampler::set_probability_and_nth_selection(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  const double nth_selection_probability = params.nth_selection != 0 ? static_cast<double>(1) / params.nth_selection : 0;
  if (nth_selection_probability > params.probability) {
    JfrSamplerWindow* const next = next_window(expired);
    // Use randomized nth selection.
    next->_sampling_interval = params.nth_selection;
    next->_nth_mod_value = randomized_nth_selection_mod_value(params.nth_selection);
    assert(next->_nth_mod_value < next->_sampling_interval, "invariant");
    return;
  }
  set_probability(params, expired);
}

void JfrAdaptiveSampler::set_all(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  const double nth_selection_probability = params.nth_selection == 0 ? 0 :
                                             static_cast<double>(1) / params.nth_selection;
  if (nth_selection_probability >= params.probability) {
    set_rate_and_nth_selection(params, expired);
    return;
  }
  set_rate_and_probability(params, expired);
}

// Called exclusively by the holder of the try lock when a window is determined to have expired.
void JfrAdaptiveSampler::rotate_window(int64_t timestamp) {
  assert(_lock, "invariant");
  EventSamplerWindow event;
  const JfrSamplerWindow* const current = active_window();
  if (!current->is_expired(timestamp)) {
    // Someone took care of it.
    return;
  }
  // debug(current, _avg_population_size);
  fill(event, current);
  rotate(current);
  event.commit();
}

inline double compute_ewma_alpha_coefficient(size_t lookback_count) {
  return lookback_count <= 1 ? 1 : static_cast<double>(1) / static_cast<double>(lookback_count);
}

inline size_t compute_accumulated_debt_carry_limit(const JfrSamplerParams& params) {
  if (params.window_duration_ms == 0 || params.window_duration_ms >= MILLIUNITS) {
    return 1;
  }
  return MILLIUNITS / params.window_duration_ms;
}

void JfrAdaptiveSampler::rotate(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  install(configure(next_window_params(expired), expired));
}

void JfrAdaptiveSampler::reset_next_window(const JfrSamplerWindow* expired) {
  next_window(expired)->reset(); // next similar window
  if (expired->is_derived()) {
    _probability_mode = false;
    next_window(expired)->reset(); // next "regular" window
    return;
  }
  _probability_mode = true;
  next_window(expired)->reset(); // next "probability" window
  _probability_mode = false;
}

const JfrSamplerWindow* JfrAdaptiveSampler::configure(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(_lock, "invariant");

  /*
   * The adaptive sampler can be configured in three dimensions (all optional):
   *
   * - rate per second  sample dynamically to maintain a continuous, maximal rate per second
   * - probability      sample using a probability
   * - nth selection    sample by selecting every nth sample point
   *
   * If a rate is specified, the probability and/or nth selection becomes relative to the rate.
   */

  static const u1 RATE = 1;
  static const u1 PROBABILITY = 2;
  static const u1 NTH_SELECTION = 4;

  #define RP  (RATE | PROBABILITY)
  #define RN  (RATE | NTH_SELECTION)
  #define PN  (PROBABILITY | NTH_SELECTION)
  #define RPN (RATE | PROBABILITY | NTH_SELECTION)

  // LEGEND
  //
  // R = "rate" option
  // P = "probability" option
  // N = "nth selection" option
  //
  // The combinations comprise an n-set (a 3-set) = { R, P, N }
  //
  // Number of r-subsets = 4 (0, 1, 2, 3) (including null set)
  //
  // Unordered selection:
  //
  // C(3, 0) = {} = NULL set = 1
  // C(3, 1) = { (R), (P), (N) } = 3
  // C(3, 2) = { (R, P), (R, N), (P, N) } = 3
  // C(3, 3) = { (R, P, N) } = 1
  //
  // in shorter terms: P({ R, P, N }) = 8
  //
  static u1 options = 0;

  if (params.reconfigure) {
    _avg_population_size = 0;
    _ewma_population_size_alpha = compute_ewma_alpha_coefficient(params.window_lookback_count);
    _acc_debt_carry_limit = compute_accumulated_debt_carry_limit(params);
    _acc_debt_carry_count = _acc_debt_carry_limit;
    reset_next_window(expired);
    assert(!_probability_mode, "invariant");
    if (params.sample_points_per_window != JfrSamplerParams::unused) {
      options |= RATE;
    }
    if (params.probability != JfrSamplerParams::unused) {
      options |= PROBABILITY;
    }
    if (params.nth_selection != JfrSamplerParams::unused) {
      options |= NTH_SELECTION;
    }
    params.reconfigure = false;
  }

  switch (options) {
    case RATE:
      set_rate(params, expired);
      break;
    case NTH_SELECTION:
      set_nth_selection(params, expired);
      break;
    case PROBABILITY:
      set_probability(params, expired);
      break;
    case RP:
      set_rate_and_probability(params, expired);
      break;
    case RN:
      set_rate_and_nth_selection(params, expired);
      break;
    case PN:
      set_probability_and_nth_selection(params, expired);
      break;
    case RPN:
      set_all(params, expired);
  }
  JfrSamplerWindow* const next = next_window(expired);
  assert(next != expired, "invariant");
  next->initialize(params, _initial_value_for_next_window);
  return next;

  #undef RP
  #undef RN
  #undef PN
  #undef RPN
}

void JfrAdaptiveSampler::install(const JfrSamplerWindow* next) {
  assert(next != NULL, "invariant");
  assert(next != active_window(), "invariant");
  Atomic::release_store(&_active_window, next);
}

inline JfrSamplerWindow* JfrAdaptiveSampler::next_window(const JfrSamplerWindow* expired) const {
  assert(expired != NULL, "invariant");
  if (expired->is_derived()) {
    if (_probability_mode) {
      // return the next probability window
      return expired == _window_2 ? _window_3 : _window_2;
    }
    // return the next "regular" window
    return expired == _window_2 ? _window_1 : _window_0;
  }
  if (_probability_mode) {
    // return the next probability window
    return expired == _window_0 ? _window_3 : _window_2;
  }
  // return the next "regular" window
  return expired == _window_0 ? _window_1 : _window_0;
}
