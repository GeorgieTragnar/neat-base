// crossover.hpp
#pragma once

#include <random>
#include <functional>
#include <memory>
#include <map>

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
    };

    explicit CrossoverOperator(const Config& config);
    core::Genome crossover(const core::Genome& parent1, const core::Genome& parent2);

private:
    using GenePtr = const std::shared_ptr<core::Gene>;
    using GenePair = std::pair<GenePtr, GenePtr>;
    using InheritanceMap = std::map<int32_t, GenePair>;

    InheritanceMap createInheritanceMap(const core::Genome& parent1, const core::Genome& parent2);
    
    Config config;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> inheritanceDist;
};

}
}
