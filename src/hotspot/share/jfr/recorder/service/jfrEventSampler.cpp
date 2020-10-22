#include "precompiled.hpp"
#include "jfr/recorder/jfrEventSetting.inline.hpp"
#include "jfr/recorder/service/jfrEventSampler.hpp"
#include "jfr/utilities/jfrAllocation.hpp"

template <typename T>
class JfrEventSamplers : public JfrCHeapObj {
 private:
  T* _samplers[LAST_EVENT_ID + 1];
 public:
  bool initialize() {
    for (int i = FIRST_EVENT_ID; i <= LAST_EVENT_ID; i++) {
      _samplers[i] = new T(static_cast<JfrEventId>(i));
      if (_samplers[i] == NULL || !_samplers[i]->initialize()) {
        return false;
      }
    }
    return true;
  }

  ~JfrEventSamplers() {
    for (int i = 0; i <= LAST_EVENT_ID; i++) {
      delete _samplers[i];
      _samplers[i] = NULL;
    }
  }

  T* get_sampler(JfrEventId event_id) {
    return _samplers[event_id];
  }
};

JfrEventSampler::JfrEventSampler(JfrEventId event_id) :
    AdaptiveSampler(80, 160),
    _event_id(event_id) {}

bool JfrEventSampler::initialize() {
  return AdaptiveSampler::initialize();
}

static JfrEventSamplers<JfrEventSampler>* _samplers = NULL;

bool JfrEventSampler::create() {
  assert(_samplers == NULL, "invariant");
  _samplers = new JfrEventSamplers<JfrEventSampler>();
  return _samplers != NULL && _samplers->initialize();
}

void JfrEventSampler::destroy() {
  if (_samplers != NULL) {
    delete _samplers;
    _samplers = NULL;
  }
}

static int64_t MIN_SAMPLES_PER_WINDOW = 20;

SamplerWindowParams JfrEventSampler::new_window_params() {
  const int64_t rate_limit = JfrEventSetting::ratelimit(_event_id);
  // a simple heuristic to derive the window size and number of samples per window from the provided rate limit
  double duration = 10; // start with 10ms window
  // duration is in milliseconds and rate_limit in samples per second
  double samples = (duration * rate_limit) / static_cast<double>(1000);
  // rebalance if too few samples are to be generated per window
  if (samples < MIN_SAMPLES_PER_WINDOW) {
    // target at least MIN_SAMPLES_PER_WINDOW samples per window
    duration *= (MIN_SAMPLES_PER_WINDOW / samples);
    samples = MIN_SAMPLES_PER_WINDOW;
  }
  SamplerWindowParams params = { static_cast<int64_t>(duration), static_cast<int64_t>(samples) };
  return params;
}

JfrEventSampler* JfrEventSampler::for_event(JfrEventId event_id) {
  assert(_samplers != NULL, "JfrEventSampler has not been properly initialized");
  return _samplers->get_sampler(event_id);
}
