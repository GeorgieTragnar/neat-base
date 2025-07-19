#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <cassert>
#include "../data/PopulationData.hpp"
#include "../data/GlobalIndexRegistry.hpp"

namespace Operator {

template<typename FitnessResultType>
std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping(const std::multimap<FitnessResultType, size_t>& fitnessResults, const std::vector<DynamicGenomeData>& genomeData, std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData, const GlobalIndexRegistry& registry);

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
        
        // Only include valid genomes (active or elite in registry and not under repair)
        auto state = registry.getState(globalIndex);
        if ((state == GenomeState::Active || state == GenomeState::Elite) && !genomeData[globalIndex].isUnderRepair) {
            // Skip species that exist in species data and are marked for elimination
            if (speciesData.find(speciesId) == speciesData.end()) {
                // Extract species ID from genome data and group indices
                speciesIndices[speciesId].push_back(globalIndex);
            }
            else if (!speciesData[speciesId].isMarkedForElimination) {
                // Extract species ID from genome data and group indices
                speciesIndices[speciesId].push_back(globalIndex);
            }
        }
        
        // Update species population size (count all genomes, valid or not)
        if (speciesData.find(speciesId) != speciesData.end()) {
            speciesData[speciesId].currentPopulationSize++;
        }
    }
    
#ifndef NDEBUG
    // Validate that all global indices in fitnessResults are unique
    std::unordered_set<size_t> uniqueGlobalIndices;
    std::vector<size_t> duplicateIndices;
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        if (uniqueGlobalIndices.find(globalIndex) != uniqueGlobalIndices.end()) {
            duplicateIndices.push_back(globalIndex);
        } else {
            uniqueGlobalIndices.insert(globalIndex);
        }
    }
    assert(duplicateIndices.empty() && "Duplicate global indices detected in fitnessResults multimap");
    
    // Validate completeness: every valid genome index should appear exactly once
    size_t totalIndices = 0;
    size_t validGenomeCount = 0;
    for (const auto& [speciesId, indices] : speciesIndices) {
        totalIndices += indices.size();
        
        // Note: Genome indices within species preserve fitness ordering from multimap,
        // not ascending index order, so no ordering constraint is needed here
    }
    
    // Count valid genomes (active or elite in registry and not under repair, and not from eliminated species)
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        auto state = registry.getState(globalIndex);
        if ((state == GenomeState::Active || state == GenomeState::Elite) && !genomeData[globalIndex].isUnderRepair) {
            // Skip species that exist in species data and are marked for elimination
            if (speciesData.find(speciesId) == speciesData.end()) {
                validGenomeCount++;
            }
            else if (!speciesData[speciesId].isMarkedForElimination) {
                validGenomeCount++;
            }
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