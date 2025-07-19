#pragma once

#include "version3/data/Genome.hpp"
#include "version3/analysis/FitnessResult.hpp"
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

namespace Operator {

class CrossoverParams;

template<typename FitnessResultType>
Genome crossover(const Genome& parentA, const FitnessResultType& fitnessA, const Genome& parentB, const FitnessResultType& fitnessB, const CrossoverParams& params);

class CrossoverParams {
public:
    CrossoverParams(double disabledGeneReactivationProbability = 0.25)
        : _disabledGeneReactivationProbability(disabledGeneReactivationProbability) {}

private:
    template<typename FitnessResultType>
    friend Genome crossover(const Genome& parentA, const FitnessResultType& fitnessA, const Genome& parentB, const FitnessResultType& fitnessB, const CrossoverParams& params);
    
    const double _disabledGeneReactivationProbability;
};

namespace {
    thread_local std::random_device rd_crossover;
    thread_local std::mt19937 gen_crossover(rd_crossover());
    thread_local std::uniform_real_distribution<double> uniform_dist_crossover(0.0, 1.0);
}

template<typename FitnessResultType>
Genome crossover(
    const Genome& parentA, 
    const FitnessResultType& fitnessA,
    const Genome& parentB, 
    const FitnessResultType& fitnessB,
    const CrossoverParams& params) {
    
    // Input validation
    const auto& nodesA = parentA.get_nodeGenes();
    const auto& nodesB = parentB.get_nodeGenes();
    const auto& connectionsA = parentA.get_connectionGenes();
    const auto& connectionsB = parentB.get_connectionGenes();
    
    // Assert I/O structure compatibility
    assert(!nodesA.empty() && !nodesB.empty());
    
    // Determine fitter parent - using direct fitness comparison
    bool parentAIsFitter = fitnessA.isBetterThan(fitnessB);
    bool equalFitness = fitnessA.isEqualTo(fitnessB);
    
    if (equalFitness) {
        // For equal fitness, randomly choose which parent contributes disjoint/excess genes
        parentAIsFitter = uniform_dist_crossover(gen_crossover) < 0.5;
    }
    
    Genome offspring = parentAIsFitter ? parentA : parentB;
    const Genome& otherParent = parentAIsFitter ? parentB : parentA;
    
    // Create lookup map for other parent's genes
    std::unordered_map<uint32_t, const NodeGene*> otherNodeMap;
    std::unordered_map<uint32_t, const ConnectionGene*> otherConnMap;
    
    for (const auto& node : otherParent.get_nodeGenes()) {
        otherNodeMap[node.get_historyID()] = &node;
    }
    for (const auto& conn : otherParent.get_connectionGenes()) {
        otherConnMap[conn.get_historyID()] = &conn;
    }

    for (auto& nodeGene : offspring.get_nodeGenes()) {
        if (otherNodeMap.count(nodeGene.get_historyID()) > 0 && uniform_dist_crossover(gen_crossover) < 0.5) {
            nodeGene.get_attributes() = otherNodeMap[nodeGene.get_historyID()]->get_attributes();
        }
    }

    for (auto& connGene : offspring.get_connectionGenes()) {
        if (otherConnMap.count(connGene.get_historyID()) > 0 && uniform_dist_crossover(gen_crossover) < 0.5) {
            connGene.get_attributes() = otherConnMap[connGene.get_historyID()]->get_attributes();
        }
        
        // Apply disabled gene reactivation probability
        if (!connGene.get_attributes().enabled && uniform_dist_crossover(gen_crossover) < params._disabledGeneReactivationProbability) {
            connGene.get_attributes().enabled = true;
        }
    }

    return offspring;
}

}