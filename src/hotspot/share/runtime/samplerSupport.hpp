#ifndef SHARE_RUNTIME_SAMPLERSUPPORT_HPP
#define SHARE_RUNTIME_SAMPLERSUPPORT_HPP

// NOTE: This code is just moved from threadHeapSampler.hpp so it can be shared between several heap samplers

class SamplerSupport : public CHeapObj<mtInternal> {
  private:

  static double fast_log2(const double& d);
  static uint64_t next_random(uint64_t rnd);

  uint64_t next_random();

  const bool _no_sync;
    // Cheap random number generator
  volatile uint64_t _rnd;

  public:
  SamplerSupport(bool no_sync = true) : _no_sync(no_sync) {
    _rnd = static_cast<uint32_t>(reinterpret_cast<uintptr_t>(this));
    if (_rnd == 0) {
      _rnd = 1;
    }
  }

  size_t next_random_geometric(size_t mean);

  double next_random_uniform();
};

#endif // SHARE_RUNTIME_SAMPLERSUPPORT_HPP
