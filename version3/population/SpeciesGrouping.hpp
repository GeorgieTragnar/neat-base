#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include "PopulationData.hpp"
#include "GlobalIndexRegistry.hpp"

namespace Population {

// Main operator function - templated on fitness result type
template<typename FitnessResultType>
std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping(
    const std::multimap<FitnessResultType, size_t>& fitnessResults,
    const std::vector<DynamicGenomeData>& genomeData,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GlobalIndexRegistry& registry
) {
    assert(!fitnessResults.empty());
    
    std::unordered_map<uint32_t, std::vector<size_t>> speciesIndices;
    
    // Reset species population sizes at start of analysis
    for (auto& [speciesId, data] : speciesData) {
        data.currentPopulationSize = 0;
    }
    
    // Single pass through fitness results: O(n) time complexity
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        
        // Only include valid genomes (active in registry and not under repair)
        if (registry.getState(globalIndex) == GenomeState::Active && !genomeData[globalIndex].isUnderRepair) {
            // Extract species ID from genome data and group indices
            speciesIndices[speciesId].push_back(globalIndex);
        }
        
        // Update species population size (count all genomes, valid or not)
        if (speciesData.find(speciesId) != speciesData.end()) {
            speciesData[speciesId].currentPopulationSize++;
        }
    }
    
#ifndef NDEBUG
    // Validate completeness: every valid genome index should appear exactly once
    size_t totalIndices = 0;
    size_t validGenomeCount = 0;
    for (const auto& [speciesId, indices] : speciesIndices) {
        totalIndices += indices.size();
        
        // Note: Genome indices within species preserve fitness ordering from multimap,
        // not ascending index order, so no ordering constraint is needed here
    }
    
    // Count valid genomes (active in registry and not under repair)
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        if (registry.getState(globalIndex) == GenomeState::Active && !genomeData[globalIndex].isUnderRepair) {
            validGenomeCount++;
        }
    }
    assert(totalIndices == validGenomeCount);
    
    // Validate no empty species vectors
    for (const auto& [speciesId, indices] : speciesIndices) {
        assert(!indices.empty());
    }
    
    // Note: New species discovered in genome data may not yet exist in speciesData
    // They will be added during dynamic data update process
#endif
    
    return speciesIndices;
}

}