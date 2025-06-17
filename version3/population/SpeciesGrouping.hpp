#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include "PopulationData.hpp"

namespace Population {

// Main operator function - templated on fitness result type
template<typename FitnessResultType>
std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping(
    const std::multimap<FitnessResultType, DynamicGenomeData>& orderedGenomeData,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
) {
    assert(!orderedGenomeData.empty());
    
    std::unordered_map<uint32_t, std::vector<size_t>> speciesIndices;
    
    // Reset species population sizes at start of analysis
    for (auto& [speciesId, data] : speciesData) {
        data.currentPopulationSize = 0;
    }
    
    // Single pass through ordered genome data: O(n) time complexity
    size_t globalIndex = 0;
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        const uint32_t speciesId = genomeData.speciesId;
        
        // Only include valid genomes (not marked for elimination or under repair)
        if (!genomeData.isMarkedForElimination && !genomeData.isUnderRepair) {
            // Extract species ID from genome data and group indices
            speciesIndices[speciesId].push_back(genomeData.genomeIndex);
        }
        
        // Update species population size (count all genomes, valid or not)
        if (speciesData.find(speciesId) != speciesData.end()) {
            speciesData[speciesId].currentPopulationSize++;
        }
        
        ++globalIndex;
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
    
    // Count valid genomes (not marked for elimination or under repair)
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        if (!genomeData.isMarkedForElimination && !genomeData.isUnderRepair) {
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