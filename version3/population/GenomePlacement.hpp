#pragma once

#include <memory>
#include <cassert>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

class GenomePlacementParams;

template<typename FitnessResultType>
void genomePlacement(
    const Genome& newGenome,
    const DynamicGenomeData& genomeData,
    uint32_t currentGeneration,
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    const GenomePlacementParams& params
);

class GenomePlacementParams {
public:
    GenomePlacementParams() = delete;
    GenomePlacementParams(size_t maxPopulationSize) 
        : _maxPopulationSize(maxPopulationSize) {
        assert(maxPopulationSize > 0 && "maxPopulationSize must be greater than 0");
    }

protected:
    template<typename FitnessResultType>
    friend void genomePlacement(
        const Genome& newGenome,
        const DynamicGenomeData& genomeData,
        uint32_t currentGeneration,
        PopulationContainer<FitnessResultType>& container,
        GlobalIndexRegistry& registry,
        const GenomePlacementParams& params
    );
    
    const size_t _maxPopulationSize;    // Maximum population size before requiring replacement
};

template<typename FitnessResultType>
inline void genomePlacement(
    const Genome& newGenome,
    const DynamicGenomeData& genomeData,
    uint32_t currentGeneration,
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    const GenomePlacementParams& params
) {
    auto logger = LOGGER("operator.GenomePlacement");
    
    // Get the current generation data
    auto& currentGenomes = container.getCurrentGenomes(currentGeneration);
    auto& currentGenomeData = container.getCurrentGenomeData(currentGeneration);
    
    // First, try to get a free index for replacement
    uint32_t targetIndex = registry.getFreeIndex();
    
    if (targetIndex == INVALID_INDEX) {
        // No free indices available - check if we can grow the population
        size_t currentSize = currentGenomes.size();
        
        if (currentSize < params._maxPopulationSize) {
            // We can grow the population
            LOG_DEBUG("genomePlacement: Adding new genome to population (size: {} -> {})", 
                     currentSize, currentSize + 1);
            targetIndex = container.push_back(currentGeneration, 
                                            Genome(newGenome), 
                                            DynamicGenomeData(genomeData));
        } else {
            // Population is at max capacity and no ReadyForReplacement indices available
            // This indicates a logic error in the state management system
            assert(false && "genomePlacement: Population at capacity but no ReadyForReplacement indices available - this indicates a state management error");
        }
    } else {
        // Using existing ReadyForReplacement slot - use strict replace() method
        LOG_DEBUG("genomePlacement: Replacing genome at index {} (reusing ReadyForReplacement slot)", targetIndex);
        
        // Verify the index is actually ready for replacement
        GenomeState state = registry.getState(targetIndex);
        assert(state == GenomeState::ReadyForReplacement && 
               "genomePlacement: Attempted to replace genome that is not marked ReadyForReplacement");
        
        container.replace(currentGeneration, targetIndex, 
                         Genome(newGenome), 
                         DynamicGenomeData(genomeData));
    }
    
    LOG_DEBUG("genomePlacement: Successfully placed genome at index {}", targetIndex);
}

}