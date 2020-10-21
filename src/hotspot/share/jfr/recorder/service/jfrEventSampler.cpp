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
      _samplers[i] = new T((JfrEventId)i);
      if (_samplers[i] == NULL) {
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

static jlong MIN_SAMPLES_PER_WINDOW = 20;
static JfrEventSamplers<JfrEventSampler>* _samplers = NULL;

bool JfrEventSampler::initialize() {
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

JfrEventSampler::JfrEventSampler(JfrEventId event_id) :
    AdaptiveSampler(80, 160),
    _event_id(event_id) {}

SamplerWindowParams JfrEventSampler::new_window_params() {
  SamplerWindowParams params;
  jlong limit = JfrEventSetting::ratelimit(_event_id);
  // a simple heuristic to derive the window size and number of samples per window from the provided rate limit
  double duration = 10; // start with 10ms window
  double samples = (duration * limit) / (double)1000; // duration is in milliseconds and limit in samples per second
  // rebalance if too few samples are to be generated per window
  if (samples < MIN_SAMPLES_PER_WINDOW) {
    // target at least MIN_SAMPLES_PER_WINDOW samples per window
    duration *= (MIN_SAMPLES_PER_WINDOW / samples);
    samples = MIN_SAMPLES_PER_WINDOW;
  }
  params.duration = (jlong)duration;
  params.sample_count = (jlong)samples;

  return params;
}

JfrEventSampler* JfrEventSampler::for_event(JfrEventId event_id) {
  assert(_samplers != NULL, "JfrEventSampler has not been properly initialized");
  return _samplers->get_sampler(event_id);
}
