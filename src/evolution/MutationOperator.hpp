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
        double weightMutationRate;
        double weightPerturbationRate;
        double newNodeRate;
        double newConnectionRate;
        double weightPerturbationRange;
        double newWeightRange;
        double activationMutationRate;
        double interTierMutationRate;
        core::ActivationGene::Config activationConfig;
        bool allowRecurrent;

        Config() = delete;
        Config(double weightMutRate,
            double weightPerturbRate,
            double nodeRate,
            double connectionRate,
            double perturbRange,
            double newWeightR,
            double actMutRate,
            double interTierRate,
            core::ActivationGene::Config actConfig,
            bool recurrent)
            : weightMutationRate(validateRate(weightMutRate, "weightMutationRate"))
            , weightPerturbationRate(validateRate(weightPerturbRate, "weightPerturbationRate"))
            , newNodeRate(validateRate(nodeRate, "newNodeRate"))
            , newConnectionRate(validateRate(connectionRate, "newConnectionRate"))
            , weightPerturbationRange(validatePositive(perturbRange, "weightPerturbationRange"))
            , newWeightRange(validatePositive(newWeightR, "newWeightRange"))
            , activationMutationRate(validateRate(actMutRate, "activationMutationRate"))
            , interTierMutationRate(validateRate(interTierRate, "interTierMutationRate"))
            , activationConfig(actConfig)
            , allowRecurrent(recurrent) {}

    private:
        static double validateRate(double value, const char* paramName) {
            if (value < 0.0 || value > 1.0) {
                throw std::invalid_argument(std::string(paramName) + " must be between 0 and 1");
            }
            return value;
        }

        static double validatePositive(double value, const char* paramName) {
            if (value <= 0.0) {
                throw std::invalid_argument(std::string(paramName) + " must be positive");
            }
            return value;
        }
    };

    explicit MutationOperator(const Config& config);
    void mutate(core::Genome& genome);

private:
    bool addNodeMutation(core::Genome& genome);
    bool addConnectionMutation(core::Genome& genome);
    void mutateActivations(core::Genome& genome);
    std::vector<std::pair<int32_t, int32_t>> findPossibleConnections(const core::Genome& genome);

    Config config;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> weightDist;    // For probability checks
    std::normal_distribution<double> perturbDist;         // For weight perturbation
    std::normal_distribution<double> newWeightDist;       // For new weight initialization
};

}
}
