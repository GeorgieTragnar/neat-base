// crossover.hpp
#pragma once

#include <random>
#include <functional>
#include <memory>
#include <map>
#include "core/ActivationGene.hpp"

namespace neat {
namespace core {
class Genome;
class Gene;
}
namespace evolution {

class CrossoverOperator {
public:
    struct Config {
        double matchingGeneInheritanceRate;
        bool inheritDisabledGenes;
        double disabledGeneReenableRate;
        double matchingActivationInheritanceRate;

        Config() = delete;
        Config(double matchingGeneRate,
            bool inheritDisabled,
            double reEnableRate,
            double matchingActivationRate)
            : matchingGeneInheritanceRate(validateRate(matchingGeneRate, "matchingGeneInheritanceRate"))
            , inheritDisabledGenes(inheritDisabled)
            , disabledGeneReenableRate(validateRate(reEnableRate, "disabledGeneReenableRate"))
            , matchingActivationInheritanceRate(validateRate(matchingActivationRate, "matchingActivationInheritanceRate"))
        {}

        private:
            static double validateRate(double value, const char* paramName) {
                if (value < 0.0 || value > 1.0) {
                    throw std::invalid_argument(std::string(paramName) + " must be between 0 and 1");
                }
                return value;
            }
    };

    explicit CrossoverOperator(const Config& config);
    core::Genome crossover(const core::Genome& parent1, const core::Genome& parent2);

private:
    using GenePtr = std::shared_ptr<core::Gene>;
    using GenePair = std::pair<GenePtr, GenePtr>;
    using InheritanceMap = std::map<int32_t, GenePair>;

    InheritanceMap createInheritanceMap(const core::Genome& parent1, const core::Genome& parent2);
    core::EActivationType inheritActivation(const core::Gene& gene1, const core::Gene& gene2, bool parent1IsFitter);
    
    Config config;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> inheritanceDist;
};

}
}
