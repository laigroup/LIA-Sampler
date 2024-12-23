#pragma once

#include <random>
#include <string>
#include <unordered_map>
#include "sampler.h"

namespace sampler {

typedef enum {
    LS,
    CDCL,
    HYBRID
}SamplingMode;

class LiaSampler : public Sampler {
    std::chrono::steady_clock::time_point time_sampling_start;
    size_t num_samples = 0;
    SamplingMode mode = LS;
    std::mt19937 mt;
    std::random_device rd;
    std::unordered_map<std::string, std::string> curr_sample;
    size_t cdcl_epoch = 1;
    double fixed_var_pct = 0.5;

    double TimeElapsed();
    void print_statistic();
   public:

    LiaSampler(z3::context* _c, std::string _smtFilePath, std::string _samplesFileDir, size_t _maxNumSamples, double _maxTimeLimit, SamplingMode _mode, unsigned seed, size_t _cdclEpoch, double _fixedVarsPct)
        : Sampler(_c, _smtFilePath, _samplesFileDir, _maxNumSamples, _maxTimeLimit), mode(_mode), cdcl_epoch(_cdclEpoch), fixed_var_pct(_fixedVarsPct){
            mt.seed(seed);
    }

    z3::tactic mk_preamble_tactic(z3::context& ctx);

    void sampling() override;
    void ls_sampling(std::ofstream& samplesFile);
    void cdcl_sampling(std::ofstream& samplesFile);
    void hybrid_sampling(std::ofstream& samplesFile);
    unsigned gen_random_seed();
};
};  // namespace sampler