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
        double matchingGeneInheritanceRate = 0.5;
        bool inheritDisabledGenes = true;
        double disabledGeneReenableRate = 0.25;
        double matchingActivationInheritanceRate = 0.5;
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
