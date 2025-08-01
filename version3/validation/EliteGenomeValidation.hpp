#pragma once

#include "version3/data/PopulationContainer.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/validation/GenomeEquals.hpp"
#include "version3/logger/Logger.hpp"
#include <cassert>

namespace Operator {

template<typename FitnessResultType>
bool validateEliteGenomes(
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& currentGeneration,
    const GlobalIndexRegistry& registry
);

// Validates that elite genomes are properly preserved between generations
template<typename FitnessResultType>
bool validateEliteGenomes(
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& currentGeneration,
    const GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("operator.EliteGenomeValidation");
    
    // Cannot validate generation 0 (no previous generation)
    if (currentGeneration == 0) {
        LOG_DEBUG("validateEliteGenomes: Skipping validation for generation 0");
        return true;
    }
    
    const uint32_t previousGeneration = currentGeneration - 1;
    
    // Get data from both generations
    const auto& prevGenomes = container.getGenomes(previousGeneration);
    const auto& prevGenomeData = container.getGenomeData(previousGeneration);
    const auto& prevFitnessResults = container.getFitnessResults(previousGeneration);
    
    const auto& currGenomes = container.getGenomes(currentGeneration);
    const auto& currGenomeData = container.getGenomeData(currentGeneration);
    const auto& currFitnessResults = container.getFitnessResults(currentGeneration);
    
    bool allValid = true;
    uint32_t eliteCount = 0;
    uint32_t validationFailures = 0;
    
    LOG_DEBUG("validateEliteGenomes: Validating elite preservation from generation {} to {}", 
             previousGeneration, currentGeneration);
    
    // Iterate through all genomes from previous generation
    for (size_t globalIndex = 0; globalIndex < prevGenomes.size(); ++globalIndex) {
        const auto& prevData = prevGenomeData[globalIndex];
        
        // Skip non-elite genomes
        if (!prevData.isElite) {
            continue;
        }
        
        eliteCount++;
        
        // Validate that global index exists in current generation
        if (globalIndex >= currGenomes.size() || globalIndex >= currGenomeData.size()) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} from generation {} not found in generation {}",
                     globalIndex, previousGeneration, currentGeneration);
            allValid = false;
            validationFailures++;
            continue;
        }
        
        const auto& prevGenome = prevGenomes[globalIndex];
        const auto& currGenome = currGenomes[globalIndex];
        const auto& currData = currGenomeData[globalIndex];
        
        // Validation 1: Genome structure must be identical
        if (!genomeEquals(prevGenome, currGenome)) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} structure differs between generations {} and {}",
                     globalIndex, previousGeneration, currentGeneration);
            allValid = false;
            validationFailures++;
        }
        
        // Validation 2: Species ID must be preserved
        if (prevData.speciesId != currData.speciesId) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} species changed from {} to {} between generations {} and {}",
                     globalIndex, prevData.speciesId, currData.speciesId, previousGeneration, currentGeneration);
            allValid = false;
            validationFailures++;
        }
        
        // Validation 3: Elite status must be preserved
        if (!currData.isElite) {
            LOG_ERROR("validateEliteGenomes: Genome {} lost elite status between generations {} and {}",
                     globalIndex, previousGeneration, currentGeneration);
            allValid = false;
            validationFailures++;
        }
        
        // Validation 4: Fitness must be identical
        // Find fitness in both generations
        FitnessResultType prevFitness{}, currFitness{};
        bool foundPrevFitness = false, foundCurrFitness = false;
        
        for (const auto& [fitness, fitnessGlobalIndex] : prevFitnessResults) {
            if (fitnessGlobalIndex == globalIndex) {
                prevFitness = fitness;
                foundPrevFitness = true;
                break;
            }
        }
        
        for (const auto& [fitness, fitnessGlobalIndex] : currFitnessResults) {
            if (fitnessGlobalIndex == globalIndex) {
                currFitness = fitness;
                foundCurrFitness = true;
                break;
            }
        }
        
        if (!foundPrevFitness) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} fitness not found in generation {}",
                     globalIndex, previousGeneration);
            allValid = false;
            validationFailures++;
        }
        
        if (!foundCurrFitness) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} fitness not found in generation {}",
                     globalIndex, currentGeneration);
            allValid = false;
            validationFailures++;
        }
        
        if (foundPrevFitness && foundCurrFitness && !prevFitness.isEqualTo(currFitness)) {
            LOG_ERROR("validateEliteGenomes: Elite genome {} fitness changed between generations {} and {}",
                     globalIndex, previousGeneration, currentGeneration);
            allValid = false;
            validationFailures++;
        }
    }
    
    if (allValid) {
        LOG_DEBUG("validateEliteGenomes: All {} elite genomes successfully validated", eliteCount);
    } else {
        LOG_ERROR("validateEliteGenomes: {} validation failures found across {} elite genomes",
                 validationFailures, eliteCount);
    }
    
    return allValid;
}

} // namespace Operator