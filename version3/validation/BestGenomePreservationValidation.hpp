#pragma once

#include "version3/data/PopulationContainer.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/population/GenerationPopulationData.hpp"
#include "version3/logger/Logger.hpp"
#include <cassert>
#include <cmath>

namespace Operator {

// Validates that the best-performing genome in each species is not marked for elimination
template<typename FitnessResultType>
void validateBestGenomePreservation(
    const GenerationPopulationData<FitnessResultType>& populationData,
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& generation,
    const GlobalIndexRegistry& registry
);

template<typename FitnessResultType>
void validateBestGenomePreservation(
    const GenerationPopulationData<FitnessResultType>& populationData,
    const PopulationContainer<FitnessResultType>& container,
    const uint32_t& generation,
    const GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("validation.BestGenomePreservation");
    
    const auto& genomeData = container.getGenomeData(generation);
    
    // Check the best genome in each species
    for (const auto& [speciesId, genomeIndices] : populationData.speciesGrouping) {
        if (genomeIndices.empty()) continue;
        
        // Get the best performer (last element in the fitness-ordered vector)
        const uint32_t bestGenomeIndex = genomeIndices.back();
        
        // Check if best genome exists in genome data
        if (bestGenomeIndex >= genomeData.size()) {
            LOG_ERROR("Best genome index {} for species {} is out of bounds", bestGenomeIndex, speciesId);
            assert(false && "Best genome index out of bounds");
        }
        
        const auto& bestGenome = genomeData[bestGenomeIndex];
        const auto genomeState = registry.getState(bestGenomeIndex);

        if (genomeState == GenomeState::ReadyForReplacement)
            continue;
        
        // Check if best genome has pending elimination counter > 0
        if (bestGenome.pendingEliminationCounter > 0) {
            LOG_ERROR("Best genome {} in species {} has pendingEliminationCounter={}", 
                     bestGenomeIndex, speciesId, bestGenome.pendingEliminationCounter);
            assert(false && "Best performing genome has elimination counter - critical algorithm error");
        }
        
        // Check if best genome is marked for elimination
        if (genomeState == GenomeState::HotElimination || genomeState == GenomeState::ColdElimination) {
            LOG_ERROR("Best genome {} in species {} is marked for elimination!", 
                     bestGenomeIndex, speciesId);
            assert(false && "Best performing genome marked for elimination - critical algorithm error");
        }
    }
    
    LOG_DEBUG("validateBestGenomePreservation: All best performers validated successfully in generation {}", 
             generation);
}

} // namespace Operator