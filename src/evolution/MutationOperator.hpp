// mutation.hpp
#pragma once
#include <random>
#include <functional>
#include <memory>
#include "core/Genome.hpp"

namespace neat {
namespace evolution {

class MutationOperator {
public:
    struct Config {
        double weightMutationRate = 0.8;
        double weightPerturbationRate = 0.9;
        double newNodeRate = 0.03;
        double newConnectionRate = 0.05;
        double weightPerturbationRange = 0.2;
        double newWeightRange = 2.0;
        bool allowRecurrent = false;
    };

    explicit MutationOperator(const Config& config);
    void mutate(core::Genome& genome);

private:
    bool addNodeMutation(core::Genome& genome);
    bool addConnectionMutation(core::Genome& genome);
    std::vector<std::pair<int32_t, int32_t>> findPossibleConnections(const core::Genome& genome);

    Config config;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> weightDist;    // For probability checks
    std::normal_distribution<double> perturbDist;         // For weight perturbation
    std::normal_distribution<double> newWeightDist;       // For new weight initialization
};

}
}
