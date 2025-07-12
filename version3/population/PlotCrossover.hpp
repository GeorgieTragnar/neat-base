#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <utility>
#include "../data/GlobalIndexRegistry.hpp"
#include "../data/PopulationData.hpp"

namespace Population {

class PlotCrossoverParams;

std::vector<std::pair<size_t, size_t>> plotCrossover(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotCrossoverParams& params,
    const GlobalIndexRegistry& registry,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
);

class PlotCrossoverParams {
public:
    PlotCrossoverParams(
        double topPerformerPercentage,
        size_t baseCrossoversPerSpecies,
        double crossoverScalingFactor,
        size_t minimumSpeciesSize
    ) : _topPerformerPercentage(topPerformerPercentage),
        _baseCrossoversPerSpecies(baseCrossoversPerSpecies),
        _crossoverScalingFactor(crossoverScalingFactor),
        _minimumSpeciesSize(minimumSpeciesSize) {
        
        assert(topPerformerPercentage > 0.0 && topPerformerPercentage <= 1.0);
        assert(baseCrossoversPerSpecies > 0);
        assert(crossoverScalingFactor >= 0.0);
        assert(minimumSpeciesSize >= 2); // Need at least 2 genomes for crossover
    }

private:
    friend std::vector<std::pair<size_t, size_t>> plotCrossover(
        const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
        const PlotCrossoverParams& params,
        const GlobalIndexRegistry& registry,
        const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
    );
    
    const double _topPerformerPercentage;
    const size_t _baseCrossoversPerSpecies;
    const double _crossoverScalingFactor;
    const size_t _minimumSpeciesSize;
};

namespace {
    thread_local std::random_device rd;
    thread_local std::mt19937 gen(rd());
}

inline std::vector<std::pair<size_t, size_t>> plotCrossover(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotCrossoverParams& params,
    const GlobalIndexRegistry& registry,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData) {
    
    assert(!speciesGroupings.empty());
    
    std::vector<std::pair<size_t, size_t>> crossoverPairs;
    
    for (const auto& [speciesId, genomeIndices] : speciesGroupings) {
        assert(!genomeIndices.empty());
        
        // Skip species marked for elimination
        auto speciesDataIt = speciesData.find(speciesId);
        if (speciesDataIt != speciesData.end() && speciesDataIt->second.isMarkedForElimination) {
            continue; // Skip this entire species
        }
        
        // Filter out genomes marked for elimination through global index registry
        std::vector<size_t> activeGenomes;
        for (size_t globalIndex : genomeIndices) {
            auto state = registry.getState(static_cast<uint32_t>(globalIndex));
            if (state == GenomeState::Active || state == GenomeState::Elite) {
                activeGenomes.push_back(globalIndex);
            }
        }
        
        // Remove duplicates from activeGenomes (in case speciesGrouping contains duplicates)
        std::sort(activeGenomes.begin(), activeGenomes.end());
        activeGenomes.erase(std::unique(activeGenomes.begin(), activeGenomes.end()), activeGenomes.end());
        
        // Skip species with no active genomes or below minimum size threshold
        size_t speciesSize = activeGenomes.size();
        if (speciesSize < params._minimumSpeciesSize) {
            continue;
        }
        
        // Calculate parent pool size from top performers (based on active genomes)
        size_t parentPoolSize = static_cast<size_t>(std::ceil(speciesSize * params._topPerformerPercentage));
        parentPoolSize = std::max(parentPoolSize, size_t(2)); // Need at least 2 parents
        parentPoolSize = std::min(parentPoolSize, speciesSize); // Can't exceed species size
        
        // Calculate number of crossovers for this species
        size_t crossoverCount = params._baseCrossoversPerSpecies + 
                               static_cast<size_t>(speciesSize * params._crossoverScalingFactor);
        
        // Generate random parent pairs
        std::uniform_int_distribution<size_t> parentDist(0, parentPoolSize - 1);
        
        for (size_t i = 0; i < crossoverCount; ++i) {
            size_t parentAIndex = parentDist(gen);
            size_t parentBIndex = parentDist(gen);
            
            // Ensure different parents (avoid self-crossover)
            while (parentBIndex == parentAIndex && parentPoolSize > 1) {
                parentBIndex = parentDist(gen);
            }
            
            // Convert to global indices (active genomes already contain global indices)
            size_t globalParentA = activeGenomes[parentAIndex];
            size_t globalParentB = activeGenomes[parentBIndex];
            
            crossoverPairs.emplace_back(globalParentA, globalParentB);
        }
    }
    
#ifndef NDEBUG
    // Validate all pairs contain valid indices
    for (const auto& [parentA, parentB] : crossoverPairs) {
        assert(parentA != parentB || crossoverPairs.size() == 1); // Allow self-crossover only if no other option
    }
#endif
    
    return crossoverPairs;
}

}