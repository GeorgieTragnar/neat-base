#pragma once

#include <vector>
#include <memory>

#include "version3/data/NodeGene.hpp"
#include "version3/data/ConnectionGene.hpp"
#include "version3/data/Genome.hpp"

namespace Operator {

class WeightMutationParams;

Genome weightMutation(const Genome& genome, const WeightMutationParams& params);

class WeightMutationParams {
public:
    enum class MutationType {
        PERTURBATION_ONLY,    // Only perturb existing weights
        REPLACEMENT_ONLY,     // Only replace weights entirely  
        MIXED                 // Both perturbation and replacement
    };
    
    WeightMutationParams() = delete;
    WeightMutationParams(double perturbationRate,
                        double replacementRate, 
                        double perturbationStrength,
                        double weightRange,
                        MutationType mutationType);

protected:
    friend Genome weightMutation(const Genome& genome, const WeightMutationParams& params);
    
    const double _perturbationRate;      // Probability to perturb a weight [0.0, 1.0]
    const double _replacementRate;       // Probability to replace a weight entirely [0.0, 1.0]
    const double _perturbationStrength;  // Standard deviation for perturbation (> 0)
    const double _weightRange;           // Range for new random weights [-range, +range] (> 0)
    const MutationType _mutationType;
};

}