#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <optional>
#include <cassert>
#include "PopulationData.hpp"

namespace Population {

// Optional parameters for optimization
struct SpeciesGroupingParams {
    // Optional capacity hints for vector pre-sizing (species_id -> expected_size)
    std::optional<std::unordered_map<uint32_t, size_t>> expectedSpeciesSizes;
    
    SpeciesGroupingParams() = default;
    
    explicit SpeciesGroupingParams(
        const std::unordered_map<uint32_t, size_t>& speciesSizes
    ) : expectedSpeciesSizes(speciesSizes) {}
};

// Core operator function - transforms ordered genome data into species-grouped indices
// Preserves the iteration order of the input map while grouping by species ID
template<typename FitnessResultType>
std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping(
    const std::map<FitnessResultType, DynamicGenomeData>& orderedGenomeData,
    const SpeciesGroupingParams& params = {}
) {
    assert(!orderedGenomeData.empty());
    
    std::unordered_map<uint32_t, std::vector<size_t>> speciesIndices;
    
    // Pre-allocate vectors if capacity hints provided
    if (params.expectedSpeciesSizes.has_value()) {
        for (const auto& [speciesId, expectedSize] : params.expectedSpeciesSizes.value()) {
            speciesIndices[speciesId].reserve(expectedSize);
        }
    }
    
    // Single pass through ordered genome data: O(n) time complexity
    size_t globalIndex = 0;
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        // Extract species ID from genome data and group indices
        speciesIndices[genomeData.speciesId].push_back(globalIndex);
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
    #endif
    
    return speciesIndices;
}

// Overload that accepts capacity hints directly
template<typename FitnessResultType>
std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping(
    const std::map<FitnessResultType, DynamicGenomeData>& orderedGenomeData,
    const std::unordered_map<uint32_t, size_t>& expectedSpeciesSizes
) {
    return speciesGrouping(orderedGenomeData, SpeciesGroupingParams(expectedSpeciesSizes));
}

// Utility function to extract species population counts from index map
std::unordered_map<uint32_t, size_t> getSpeciesPopulationCounts(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesIndices
) {
    std::unordered_map<uint32_t, size_t> populationCounts;
    populationCounts.reserve(speciesIndices.size());
    
    for (const auto& [speciesId, indices] : speciesIndices) {
        populationCounts[speciesId] = indices.size();
    }
    
    return populationCounts;
}

// Utility function to resolve global index from species relative index
size_t resolveGlobalIndex(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesIndices,
    uint32_t speciesId,
    size_t relativeIndex
) {
    assert(speciesIndices.find(speciesId) != speciesIndices.end());
    assert(relativeIndex < speciesIndices.at(speciesId).size());
    
    return speciesIndices.at(speciesId)[relativeIndex];
}

// Utility function to get the total number of genomes across all species
size_t getTotalGenomeCount(const std::unordered_map<uint32_t, std::vector<size_t>>& speciesIndices) {
    size_t totalCount = 0;
    for (const auto& [speciesId, indices] : speciesIndices) {
        totalCount += indices.size();
    }
    return totalCount;
}

}