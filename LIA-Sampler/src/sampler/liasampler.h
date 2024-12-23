#pragma once

#include "sampler.h"

namespace sampler {
class LiaSampler : public Sampler {
    std::chrono::steady_clock::time_point time_sampling_start;
    size_t num_samples = 0;

    double TimeElapsed();
    void print_statistic();
   public:
    LiaSampler(z3::context* _c, std::string _smtFilePath, std::string _samplesFileDir, size_t _maxNumSamples, double _maxTimeLimit)
        : Sampler(_c, _smtFilePath, _samplesFileDir, _maxNumSamples, _maxTimeLimit) {
    }

    z3::tactic mk_preamble_tactic(z3::context& ctx);

    void sampling() override;
};
};  // namespace sampler