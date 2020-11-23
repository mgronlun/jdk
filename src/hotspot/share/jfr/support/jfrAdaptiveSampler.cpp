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
  return Thread::current()->jfr_thread_local()->next_random_uniform();
}

class JfrProbabilitySamplerWindow : public JfrSamplerWindow {
  friend class JfrAdaptiveSampler;
 private:
  double _probability;

  JfrProbabilitySamplerWindow() : JfrSamplerWindow(), _probability(MAX_PROBABILITY) {}

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
   * For continuous probabilities, the decision to sample
   * is a function of a Bernoulli trial, where the continuous
   * random variable X takes a uniformly distributed value
   * in the interval [0,1]. The probability value assigned
   * to the window is the p value of the trial, i.e. the
   * probability of success. If a rate is specified, then
   * the attempt must also be below the sample size maximum.
   *
   * If X < p, the trial outcome is successful.
   */
  bool sample() const override {
    const double p = probability();
    if (p == MAX_PROBABILITY) return true;
    const size_t max = max_sample_size();
    if (population_size() > max || next_random_uniform() >= p) {
      // Event did not occur.
      return false;
    }
    return Atomic::add(&_measured_population_size, static_cast<size_t>(1)) <= max;
  }

 public:
  size_t sample_size() const override {
    const size_t size = population_size();
    const size_t max_size = max_sample_size();
    return size < max_size ? size : max_size;
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
  _measured_population_size(0) {
  _params.probability = JfrSamplerParams::unused;
  _params.nth_selection = JfrSamplerParams::unused;
}

void JfrSamplerWindow::initialize(const JfrSamplerParams& params, size_t initializer) {
  _params = params;
  if (params.window_duration_ms == 0) {
    Atomic::store(&_end_ticks, static_cast<int64_t>(0));
    return;
  }
  Atomic::store(&_measured_population_size, initializer);
  const int64_t end_ticks = now() + millis_to_countertime(params.window_duration_ms);
  Atomic::store(&_end_ticks, end_ticks);
}

void JfrSamplerWindow::reset() {
  _nth_mod_value = 0;
  _sampling_interval = 1;
  _projected_sample_size = MAX_SIZE_T;
  _projected_population_size = MAX_SIZE_T;
}

inline size_t JfrSamplerWindow::population_size() const {
  return Atomic::load(&_measured_population_size);
}

inline size_t JfrSamplerWindow::max_sample_size() const {
  return _projected_sample_size;
}

// For nth selections, representing discrete probabilities,
// the sample size is derived from the measured population size.
size_t JfrSamplerWindow::sample_size() const {
  const size_t size = population_size();
  if (size > _projected_population_size) {
    return max_sample_size();
  }
  if (size == 0 || size < _nth_mod_value) {
    return 0;
  }
  if (_nth_mod_value == 0) {
    return size / _sampling_interval;
  }
  return ((size - _nth_mod_value) / _sampling_interval) + 1;
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
  _next_window_initializer(0),
  _lock(0) {}

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
  _window_2 = new JfrProbabilitySamplerWindow();
  if (_window_2 == NULL) {
    return false;
  }
  assert(_window_3 == NULL, "invariant");
  _window_3 = new JfrProbabilitySamplerWindow();
  if (_window_3 == NULL) {
    return false;
  }
  _active_window = _window_0;
  return true;
}

inline bool JfrSamplerWindow::is_derived() const {
  return false;
}

// Maps an expired window to the next window to be installed, as a function of
// the type of the just expired window and the configuration of the sampler.
inline JfrSamplerWindow* JfrAdaptiveSampler::next_window(const JfrSamplerWindow* expired, bool probability /* false */) const {
  assert(expired != NULL, "invariant");
  if (expired->is_derived()) {
    if (probability) {
      // return the next probability window
      return expired == _window_2 ? _window_3 : _window_2;
    }
    // return the next "regular" window
    return expired == _window_2 ? _window_1 : _window_0;
  }
  if (probability) {
    // return the next probability window
    return expired == _window_0 ? _window_3 : _window_2;
  }
  // return the next "regular" window
  return expired == _window_0 ? _window_1 : _window_0;
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

// Randomly select the nth element in the interval, map it to its zero based index.
inline size_t randomized_nth_selection_mod_value(size_t interval) {
  assert(interval > 0, "invariant");
  if (interval <= 1) {
    return 0;
  }
  const double factor = next_random_uniform();
  assert(factor >= 0.0, "invariant");
  assert(factor <= 1.0, "invariant");
  return round(factor * (interval - 1));
}

/*
 * Based on previous information, the sampler creates a future 'projection',
 * a speculation of what the situation will be like for the next window.
 * Given this projection, the parameters are set accordingly to collect
 * a sample set that is as close as possible to the target (set point), which
 * is a function of the number of sample_points_per_window + accumulated debt.
 */
JfrSamplerWindow* JfrAdaptiveSampler::set_rate(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.probability == JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection == JfrSamplerParams::unused, "invariant");
  JfrSamplerWindow* const next = set_projected_sample_size(params, expired);
  if (next->_projected_sample_size == 0) {
    return next;
  }
  // _avg_population_size = next_window_population_size(expired, _ewma_population_size_alpha, _avg_population_size);
  if (expired->population_size() > next->_projected_sample_size) {
    next->_sampling_interval = expired->population_size() / next->_projected_sample_size;
    // next->_nth_mod_value = randomized_nth_selection_mod_value(next->_sampling_interval);
  } else {
    next->_sampling_interval = 1;
    //next->_nth_mod_value = 0;
  }
  next->_nth_mod_value = 0;
  return set_projected_population_size(next->_projected_sample_size, next);
  // next->_projected_population_size = next->_projected_sample_size * next->_sampling_interval;
  //return next;
}

JfrSamplerWindow* JfrAdaptiveSampler::set_projected_sample_size(const JfrSamplerParams& params, const JfrSamplerWindow* expired, bool probability /* false */) {
  assert(params.sample_points_per_window != JfrSamplerParams::JfrSamplerParams::unused, "invariant");
  JfrSamplerWindow* const next = next_window(expired, probability);
  assert(next != NULL, "invariant");
  assert(next != expired, "invariant");
  next->_projected_sample_size = params.sample_points_per_window + next_window_amortization(expired);
  return next;
}

/*
 * When the sampler is configured to maintain a rate, is employs the concepts
 * of 'debt' and 'accumulated debt'. 'Accumulated debt' can be thought of as
 * a cumulative error term, and is indicative for how much the sampler is
 * deviating from a set point, i.e. the target rate. Debt accumulates naturally
 * as a function of undersampled windows, mainly because of system fluctuations,
 * i.e. too small populations.
 *
 * A specified rate is implicitly a _maximal_ rate, so the sampler must ensure
 * to respect this 'limit'. Rates are normalized as per-second ratios, hence the
 * limit to respect is on a per second basis. During this second, the sampler
 * has freedom to dynamically re-adjust, and it does so by 'amortizing'
 * accumulated debt over a certain number of windows that fall within the second.
 *
 * Intuitively, accumulated debt 'carry over' from the predecessor to the successor
 * window if within the allowable time frame (determined in # of 'windows' given by
 * _acc_debt_carry_limit). The successor window will sample more points to make amends,
 * or 'amortize' debt accumulated by its predecessor(s).
 */
size_t JfrAdaptiveSampler::next_window_amortization(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  return -expired->accumulated_debt();
  /*
  if (_acc_debt_carry_count == _acc_debt_carry_limit) {
    _acc_debt_carry_count = 1;
    return 0;
  }
  ++_acc_debt_carry_count;
  return -expired->accumulated_debt(); // negation
  */
}

inline intptr_t JfrSamplerWindow::debt() const {
  return static_cast<intptr_t>(sample_size() - _params.sample_points_per_window);
}

inline intptr_t JfrSamplerWindow::accumulated_debt() const {
  return static_cast<intptr_t>(_params.sample_points_per_window - max_sample_size()) + debt();
}

JfrSamplerWindow* JfrAdaptiveSampler::set_projected_population_size(size_t projected_sample_size, JfrSamplerWindow* next) {
  assert(projected_sample_size > 0, "invariant");
  assert(next != NULL, "invariant");
  assert(next != active_window(), "invariant");
  assert(next->_sampling_interval > 0, "invariant");
  assert(next->_nth_mod_value < next->_sampling_interval, "invariant");
  projected_sample_size -= next->_nth_mod_value != 0 ? 1 : 0;
  next->_projected_population_size = (projected_sample_size * next->_sampling_interval) + next->_nth_mod_value;
  return next;
}

JfrSamplerWindow* JfrAdaptiveSampler::set_probability(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  assert(params.sample_points_per_window == JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection == JfrSamplerParams::unused, "invariant");
  JfrProbabilitySamplerWindow* const next = static_cast<JfrProbabilitySamplerWindow*>(next_window(expired, true));
  assert(next != expired, "invariant");
  if (params.probability == 0.0) {
    next->_projected_sample_size = 0;
    return next;
  }
  next->_projected_sample_size = MAX_SIZE_T;
  next->set_probability(params.probability);
  return next;
}

JfrSamplerWindow* JfrAdaptiveSampler::set_rate_and_probability(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection == JfrSamplerParams::unused, "invariant");
  JfrProbabilitySamplerWindow* const next = static_cast<JfrProbabilitySamplerWindow*>(set_projected_sample_size(params, expired, true));
  assert(next != expired, "invariant");
  if (next->_projected_sample_size == 0) {
    return next;
  }
  next->set_probability(params.probability);
  return next;
}

JfrSamplerWindow* JfrAdaptiveSampler::set_nth_selection(const JfrSamplerParams& params, const JfrSamplerWindow* expired, bool randomize) {
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  assert(params.sample_points_per_window == JfrSamplerParams::unused, "invariant");
  assert(params.probability == JfrSamplerParams::unused, "invariant");
  // To keep the modulo consistent across window boundaries.
  _next_window_initializer = expired->population_size() % expired->_sampling_interval;
  JfrSamplerWindow* const next = next_window(expired);
  assert(next != expired, "invariant");
  next->_projected_population_size = MAX_SIZE_T;
  next->_sampling_interval = params.nth_selection;
  next->_nth_mod_value = randomize ? randomized_nth_selection_mod_value(params.nth_selection) : 0;
  return next;
}

JfrSamplerWindow* JfrAdaptiveSampler::set_rate_and_nth_selection(const JfrSamplerParams& params, const JfrSamplerWindow* expired, bool randomize) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
  assert(params.probability == JfrSamplerParams::unused, "invariant");
  JfrSamplerWindow* const next = set_projected_sample_size(params, expired);
  if (next->_projected_sample_size == 0) {
    return next;
  }
  next->_sampling_interval = params.nth_selection;
  next->_nth_mod_value = randomize ? randomized_nth_selection_mod_value(params.nth_selection) : 0;
  return set_projected_population_size(next->_projected_sample_size, next);
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

void JfrAdaptiveSampler::install(const JfrSamplerWindow* next) {
  assert(next != NULL, "invariant");
  assert(next != active_window(), "invariant");
  Atomic::release_store(&_active_window, next);
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

void JfrAdaptiveSampler::rotate(const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(expired == active_window(), "invariant");
  install(configure(next_window_params(expired), expired));
}

static constexpr const u1 RATE = 1;
static constexpr const u1 PROBABILITY = 2;
static constexpr const u1 NTH_SELECTION = 4;
static constexpr const u1 RATE_PROBABIILTY = RATE | PROBABILITY;
static constexpr const u1 RATE_NTH_SELECTION = RATE | NTH_SELECTION;

const JfrSamplerWindow* JfrAdaptiveSampler::configure(const JfrSamplerParams& params, const JfrSamplerWindow* expired) {
  assert(expired != NULL, "invariant");
  assert(_lock, "invariant");
  //
  // The adaptive sampler can be configured in three dimensions (all optional):
  //
  // - rate per second  sample dynamically to maintain a continuous, maximal rate per second
  // - probability      sample using a probability
  // - nth selection    sample by selecting every nth sample point
  //
  //
  // LEGEND
  //
  // R = "rate" option
  // P = "probability" option
  // N = "nth selection" option
  //
  // The combinations comprise an n-set (a 3-set) = { R, P, N }, but the options for probability
  // and nth selection are made mutually exclusive or disjoint after preprocessing.
  // Therefore, the semantically valid combinations constitute only a strict subset of the power set.
  //
  // Number of r-subsets = 3 (0, 1, 2) (including null set)
  //
  // Unordered selection:
  //
  // C(3, 0) = {} = NULL set = 1
  // C(3, 1) = { (R), (P), (N) } = 3
  // C(3, 2) = { (R, P), (R, N) } = 2
  //
  // | P({ R, P, N }) \ { (P, N), (R, P, N) } | = 6
  //
  static u1 subset = 0;
  static bool randomize_nth_selection = false;
  if (params.reconfigure) {
    subset = configure(params, expired, &randomize_nth_selection);
    assert(!params.reconfigure, "invariant");
  }
  JfrSamplerWindow* next = NULL;
  switch (subset) {
    case RATE:
      next = set_rate(params, expired);
      break;
    case PROBABILITY:
      next = set_probability(params, expired);
      break;
    case NTH_SELECTION:
      next = set_nth_selection(params, expired, randomize_nth_selection);
      break;
    case RATE_PROBABIILTY:
      next = set_rate_and_probability(params, expired);
      break;
    case RATE_NTH_SELECTION:
      next = set_rate_and_nth_selection(params, expired, randomize_nth_selection);
      break;
    case 0: // NULL set
    default:
      next = next_window(expired);
  }
  assert(next != NULL, "invariant");
  assert(next != expired, "invariant");
  next->initialize(params, _next_window_initializer);
  return next;
}

void JfrAdaptiveSampler::reset_next_window(const JfrSamplerWindow* expired) {
  next_window(expired)->reset(); // next isomorphic window
  if (expired->is_derived()) {
    next_window(expired)->reset(); // next "regular" window
    return;
  }
  next_window(expired, true)->reset(); // next "probability" window
}

u1 JfrAdaptiveSampler::configure(const JfrSamplerParams& params, const JfrSamplerWindow* expired, bool* randomize_nth_selection) {
  assert(params.reconfigure, "invariant");
  assert(expired != NULL, "invariant");
  assert(randomize_nth_selection != NULL, "invariant");
  _next_window_initializer = 0;
  reset_next_window(expired);
  u1 subset = 0;
  *randomize_nth_selection = false;
  if (params.sample_points_per_window != JfrSamplerParams::unused) {
    subset |= configure_rate(params);
  }
  if (params.probability != JfrSamplerParams::unused) {
    subset |= configure_probability(params, randomize_nth_selection);
  }
  if (params.nth_selection != JfrSamplerParams::unused) {
    assert(params.probability == JfrSamplerParams::unused, "not mutually exclusive");
    subset |= NTH_SELECTION;
  }
  params.reconfigure = false;
  return subset;
}

u1 JfrAdaptiveSampler::configure_rate(const JfrSamplerParams& params) {
  assert(params.sample_points_per_window != JfrSamplerParams::unused, "invariant");
  _avg_population_size = 0;
  _ewma_population_size_alpha = compute_ewma_alpha_coefficient(params.window_lookback_count);
  _acc_debt_carry_limit = compute_accumulated_debt_carry_limit(params);
  _acc_debt_carry_count = _acc_debt_carry_limit;
  return RATE;
}

u1 JfrAdaptiveSampler::configure_probability(const JfrSamplerParams& params, bool* randomize_nth_selection) const {
  if (transform_probability_to_nth_selection(const_cast<JfrSamplerParams&>(params), randomize_nth_selection)) {
    assert(params.probability == JfrSamplerParams::unused, "invariant");
    assert(params.nth_selection != JfrSamplerParams::unused, "invariant");
    // we'll use nth selection
    return 0;
  }
  assert(params.nth_selection == JfrSamplerParams::unused, "invariant");
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  return PROBABILITY;
}

static bool is_uniform_sample_space_discrete(double p, size_t* nth_selection) {
  assert(nth_selection != NULL, "invariant");
  assert(p >= 0.0, "invariant");
  assert(p <= MAX_PROBABILITY, "invariant");
  if (p == 0.0) {
    *nth_selection = 0;
    return true;
  }
  const double sample_space_size = static_cast<double>(1) / p;
  *nth_selection = sample_space_size;
  return sample_space_size == static_cast<double>(*nth_selection);
}

bool JfrAdaptiveSampler::transform_probability_to_nth_selection(JfrSamplerParams& params, bool* randomize_nth_selection) const {
  assert(params.probability != JfrSamplerParams::unused, "invariant");
  assert(randomize_nth_selection != NULL, "invariant");
  assert(!*randomize_nth_selection, "invariant");
  if (params.nth_selection != JfrSamplerParams::unused) {
    const double nth_selection_probability = params.nth_selection == 0 ? 0 : static_cast<double>(1) / params.nth_selection;
    if (nth_selection_probability > params.probability) {
      // covered by nth selection
      params.probability = JfrSamplerParams::unused;
      *randomize_nth_selection = false;
      return true;
    }
    // nth selection covered by probability
    params.nth_selection = JfrSamplerParams::unused;
  }
  size_t probability_nth_selection;
  if (is_uniform_sample_space_discrete(params.probability, &probability_nth_selection)) {
    params.nth_selection = static_cast<intptr_t>(probability_nth_selection);
    params.probability = JfrSamplerParams::unused;
    *randomize_nth_selection = true;
    return true;
  }
  return false;
}
