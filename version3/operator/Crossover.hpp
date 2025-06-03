#pragma once

#include "version3/data/Genome.hpp"
#include "version3/analysis/FitnessResult.hpp"
#include "version3/analysis/GenomeAnalytics.hpp"

namespace Operator {

class CrossoverParams;

Genome crossover(const Genome& parentA, const Analysis::GenomeAnalytics& analyticsA, const Genome& parentB, const Analysis::GenomeAnalytics& analyticsB, const CrossoverParams& params);

class CrossoverParams {
public:
    CrossoverParams(double disabledGeneReactivationProbability = 0.25)
        : _disabledGeneReactivationProbability(disabledGeneReactivationProbability) {}

private:
    friend Genome crossover(const Genome& parentA, const Analysis::GenomeAnalytics& analyticsA, const Genome& parentB, const Analysis::GenomeAnalytics& analyticsB, const CrossoverParams& params);
    
    const double _disabledGeneReactivationProbability;
};

}