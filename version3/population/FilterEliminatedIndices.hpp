#pragma once

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>
#include "../data/GlobalIndexRegistry.hpp"
#include "GenerationPopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

// FilterEliminatedIndices operator: Filters species grouping in-place to exclude eliminated genomes
template<typename FitnessResultType>
void filterEliminatedIndices(
    GenerationPopulationData<FitnessResultType>& populationData,
    const GlobalIndexRegistry& registry
);

// Template implementation
template<typename FitnessResultType>
void filterEliminatedIndices(
    GenerationPopulationData<FitnessResultType>& populationData,
    const GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("population.FilterEliminatedIndices");
    
    size_t originalTotalGenomes = 0;
    size_t originalSpeciesCount = populationData.speciesGrouping.size();
    
    // Count original genomes for logging
    for (const auto& [speciesId, genomeIndices] : populationData.speciesGrouping) {
        originalTotalGenomes += genomeIndices.size();
    }
    
    // Phase 1: Filter eliminated genomes from each species in-place
    for (auto& [speciesId, genomeIndices] : populationData.speciesGrouping) {
        // Remove eliminated genomes in-place using std::remove_if
        auto newEnd = std::remove_if(genomeIndices.begin(), genomeIndices.end(),
            [&registry](uint32_t globalIndex) {
                return registry.getState(globalIndex) != GenomeState::Active;
            });
        
        // Erase the eliminated elements
        genomeIndices.erase(newEnd, genomeIndices.end());
    }
    
    // Phase 2: Remove species that have no active genomes
    for (auto it = populationData.speciesGrouping.begin(); 
         it != populationData.speciesGrouping.end();) {
        if (it->second.empty()) {
            LOG_DEBUG("Species {} has no active genomes after filtering - removed from grouping", it->first);
            it = populationData.speciesGrouping.erase(it);
        } else {
            ++it;
        }
    }
    
    // Count final genomes for logging
    size_t finalTotalGenomes = 0;
    for (const auto& [speciesId, genomeIndices] : populationData.speciesGrouping) {
        finalTotalGenomes += genomeIndices.size();
    }
    
    LOG_DEBUG("Filtered species grouping in-place: {} → {} genomes, {} → {} species",
             originalTotalGenomes, finalTotalGenomes, originalSpeciesCount, populationData.speciesGrouping.size());
    
#ifndef NDEBUG
    // Validation: ensure no empty species vectors remain
    for (const auto& [speciesId, indices] : populationData.speciesGrouping) {
        assert(!indices.empty() && "Filtered species grouping should not contain empty species vectors");
        
        // Validation: ensure all indices are actually active
        for (uint32_t globalIndex : indices) {
            assert(registry.getState(globalIndex) == GenomeState::Active && 
                   "Filtered grouping should only contain active genomes");
        }
    }
#endif
}

} // namespace Operator