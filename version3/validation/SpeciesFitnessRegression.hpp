#pragma once

#include "version3/data/PopulationContainer.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/logger/Logger.hpp"
#include <unordered_map>
#include <cassert>

namespace Operator {

#ifndef NDEBUG

// Forward declaration
template<typename FitnessResultType>
void validateSpeciesFitnessRegression(
    std::unordered_map<uint32_t, FitnessResultType>& speciesBestFitnessTracker,
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& generation,
    const GlobalIndexRegistry& registry
);

// Validates that species fitness doesn't regress and updates the tracking map
template<typename FitnessResultType>
void validateSpeciesFitnessRegression(
    std::unordered_map<uint32_t, FitnessResultType>& speciesBestFitnessTracker,
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& generation,
    const GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("operator.SpeciesFitnessRegression");
    
    // Get current generation data
    const auto& genomeData = container.getGenomeData(generation);
    const auto& fitnessResults = container.getFitnessResults(generation);
    
    // Step 1: Build temporary map of current generation's best fitness per species
    std::unordered_map<uint32_t, FitnessResultType> currentSpeciesBest;
    
    // Iterate through fitness results (ordered from worst to best)
    // The last (best) fitness seen for each species will remain in the map
    for (const auto& [fitness, globalIndex] : fitnessResults) {
        // Skip if genome is not Active
        if (registry.getState(globalIndex) != GenomeState::Active) {
            continue;
        }
        
        // Get species ID for this genome
        if (globalIndex < genomeData.size()) {
            uint32_t speciesId = genomeData[globalIndex].speciesId;
            
            // Update with better fitness (overwrites previous entry)
            currentSpeciesBest[speciesId] = fitness;
        }
    }
    
    LOG_DEBUG("validateSpeciesFitnessRegression: Generation {} - Tracking {} species", 
             generation, currentSpeciesBest.size());
    
    // Step 2: Validate against historical best and update tracker
    for (const auto& speciesBest : currentSpeciesBest) {
        const uint32_t speciesId = speciesBest.first;
        const FitnessResultType& currentBest = speciesBest.second;
        
        auto trackerIt = speciesBestFitnessTracker.find(speciesId);
        
        if (trackerIt != speciesBestFitnessTracker.end()) {
            // Species exists in tracker - validate no regression
            const FitnessResultType& historicalBest = trackerIt->second;
            
            // Assert that current best is not worse than historical best
            if (historicalBest.isBetterThan(currentBest)) {
                LOG_ERROR("Species {} fitness regression detected in generation {}: current best < historical best", 
                         speciesId, generation);
                assert(false && "Species fitness regression detected - evolution should not lose fitness gains");
            }
            
            // Update tracker if current is better
            if (currentBest.isBetterThan(historicalBest)) {
                LOG_DEBUG("Species {} improved fitness in generation {} (updating tracker)", 
                         speciesId, generation);
                trackerIt->second = currentBest;
            }
        } else {
            // New species - record its first best fitness
            LOG_DEBUG("Species {} discovered in generation {} (adding to tracker)", 
                     speciesId, generation);
            speciesBestFitnessTracker[speciesId] = currentBest;
        }
    }
    
    LOG_DEBUG("validateSpeciesFitnessRegression: Tracker now contains {} species", 
             speciesBestFitnessTracker.size());
}

#endif // NDEBUG

} // namespace Operator