#pragma once

#include <memory>
#include <cassert>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/PopulationData.hpp"
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
    GenomePlacementParams(
        size_t maxPopulationSize,
        bool protectElites = true
    ) : _maxPopulationSize(maxPopulationSize),
        _protectElites(protectElites) {
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
    const bool _protectElites;          // Whether to protect elite genomes from replacement
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
        // No free indices available - check if we need to expand or if we're at capacity
        size_t currentSize = currentGenomes.size();
        
        if (currentSize < params._maxPopulationSize) {
            // We can grow the population
            LOG_DEBUG("genomePlacement: Adding new genome to population (size: {} -> {})", 
                     currentSize, currentSize + 1);
            targetIndex = container.push_back(currentGeneration, 
                                            Genome(newGenome), 
                                            DynamicGenomeData(genomeData));
        } else {
            // Population is at max capacity - we need to force replacement
            // TODO: Implement smarter replacement strategy (e.g., replace worst fitness, oldest genome, etc.)
            // TODO: Consider replacement policies that preserve species diversity
            // TODO: Add metrics/logging for forced replacement events
            
            // Find a non-elite genome to replace
            uint32_t replaceIndex = INVALID_INDEX;
            
            for (uint32_t i = 0; i < currentSize; ++i) {
                GenomeState state = registry.getState(i);
                if (!params._protectElites || state != GenomeState::Elite) {
                    // TODO: Instead of taking first non-elite, implement selection criteria
                    // (e.g., worst fitness, oldest age, furthest from species representative)
                    replaceIndex = i;
                    break;
                }
            }
            
            if (replaceIndex == INVALID_INDEX) {
                // All genomes are elite and protected - cannot place new genome
                // TODO: Consider policy options: force replace anyway, queue for later, or reject
                // TODO: Add proper error handling/return mechanism to signal placement failure
                LOG_WARN("genomePlacement: Cannot place genome - population at capacity and all genomes are protected elites");
                return;
            }
            
            LOG_DEBUG("genomePlacement: Force replacing genome {} (population at capacity: {})", 
                     replaceIndex, currentSize);
            
            // Force replace the selected genome
            currentGenomes[replaceIndex] = Genome(newGenome);
            currentGenomeData[replaceIndex] = genomeData;
            
            // TODO: Proper state management for forced replacement
            // The state might need to be reset more carefully depending on what was replaced
            // TODO: Consider updating elimination counters or other metadata when forcing replacement
            // TODO: Handle potential race conditions if multiple threads are placing genomes
            
            targetIndex = replaceIndex;
        }
    } else {
        // Using existing eliminated slot - use strict replace() method
        LOG_DEBUG("genomePlacement: Replacing genome at index {} (reusing eliminated slot)", targetIndex);
        container.replace(currentGeneration, targetIndex, 
                         Genome(newGenome), 
                         DynamicGenomeData(genomeData));
    }
    
    // TODO: Add validation that placement was successful and data is consistent
    // TODO: Consider adding placement statistics/metrics for monitoring
    // TODO: Handle fitness results - should this operator also manage fitness data placement?
    
    LOG_DEBUG("genomePlacement: Successfully placed genome at index {}", targetIndex);
}

}