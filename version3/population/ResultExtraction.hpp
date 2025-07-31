#pragma once

#include <cassert>
#include "../data/PopulationContainer.hpp"
#include "../data/GlobalIndexRegistry.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

template<typename FitnessResultType>
uint32_t extractBestGenomeIndex(
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const GlobalIndexRegistry& registry
);

template<typename FitnessResultType>
uint32_t extractBestGenomeIndex(
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("population.ResultExtraction");
    
    const auto& fitnessResults = container.getFitnessResults(generation);
    
    assert(!fitnessResults.empty() && "extractBestGenomeIndex: No fitness results found");
    
    LOG_DEBUG("extractBestGenomeIndex: Searching for best Active genome from {} fitness results", fitnessResults.size());
    
    // Search backwards through fitness results (highest fitness first) for best Active genome
    for (auto it = fitnessResults.rbegin(); it != fitnessResults.rend(); ++it) {
        uint32_t candidateIndex = it->second;
        
        // Check if this genome is in Active state
        if (registry.getState(candidateIndex) == GenomeState::Active) {
            LOG_DEBUG("extractBestGenomeIndex: Found best Active genome at index {}", candidateIndex);
            return candidateIndex;
        }
    }
    
    // If we reach here, no Active genomes were found - this is a critical error
    assert(false && "extractBestGenomeIndex: No Active genomes found in population");
    return UINT32_MAX; // Should never reach here due to assertion
}

} // namespace Operator