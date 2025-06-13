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
        
        // Extract species ID from genome data and group indices
        speciesIndices[speciesId].push_back(globalIndex);
        
        // Update species population size
        if (speciesData.find(speciesId) != speciesData.end()) {
            speciesData[speciesId].currentPopulationSize++;
        }
        
        ++globalIndex;
    }
    
    #ifdef DEBUG
    // Validate completeness: every global index should appear exactly once
    size_t totalIndices = 0;
    for (const auto& [speciesId, indices] : speciesIndices) {
        totalIndices += indices.size();
        
        // Validate ordering preservation within species
        for (size_t i = 1; i < indices.size(); ++i) {
            assert(indices[i-1] < indices[i]);
        }
    }
    assert(totalIndices == orderedGenomeData.size());
    
    // Validate no empty species vectors
    for (const auto& [speciesId, indices] : speciesIndices) {
        assert(!indices.empty());
    }
    
    // Validate all species in result exist in species data
    for (const auto& [speciesId, indices] : speciesIndices) {
        assert(speciesData.find(speciesId) != speciesData.end());
    }
    #endif
    
    return speciesIndices;
}

}