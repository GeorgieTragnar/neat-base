#pragma once

#include <memory>
#include <cassert>
#include <functional>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/analysis/FitnessStrategy.hpp"
#include "version3/analysis/SpeciationControlUnit.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

template<typename FitnessResultType>
void genomePlacement(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    Genome genome,
    std::function<DynamicGenomeData(const Genome&, uint32_t)> metadataCreator,
    uint32_t currentGeneration,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl
);

template<typename FitnessResultType>
inline void genomePlacement(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    Genome genome,
    std::function<DynamicGenomeData(const Genome&, uint32_t)> metadataCreator,
    uint32_t currentGeneration,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl
) {
    auto logger = LOGGER("operator.GenomePlacement");
    
    // First, try to get a free index for replacement
    uint32_t targetIndex = registry.getFreeIndex();
    
    if (targetIndex == INVALID_INDEX) {
        // No free indices available - add new genome to population
        // CRITICAL FIX: Use the REAL global index returned by push_back, not vector size
        DynamicGenomeData placeholderMetadata{};
        targetIndex = container.push_back(currentGeneration, std::move(genome), std::move(placeholderMetadata));
        LOG_DEBUG("genomePlacement: Adding new genome to population at global index {}", targetIndex);
        
        // Get reference to placed genome and create actual metadata
        auto& genomes = container.getGenomes(currentGeneration);
        auto& genomeData = container.getGenomeData(currentGeneration);
        DynamicGenomeData metadata = metadataCreator(genomes[targetIndex], targetIndex);
        
        // Update with actual metadata
        genomeData[targetIndex] = std::move(metadata);
        
        // Perform fitness evaluation if genome is evaluable
        if (!genomeData[targetIndex].isUnderRepair) {
            FitnessResultType fitness = fitnessStrategy->evaluate(
                genomes[targetIndex].get_phenotype(), 
                speciationControl
            );
            
            // Insert fitness result into multimap
            auto& fitnessResults = container.getFitnessResults(currentGeneration);
            fitnessResults.insert({fitness, targetIndex});
        }
        
    } else {
        LOG_DEBUG("genomePlacement: Replacing genome at global index {} (reusing ReadyForReplacement slot)", targetIndex);
        
        // Verify the index is actually ready for replacement
        GenomeState state = registry.getState(targetIndex);
        
        // Replace genome in container first with placeholder metadata
        auto& genomes = container.getGenomes(currentGeneration);
        auto& genomeData = container.getGenomeData(currentGeneration);
        genomes[targetIndex] = std::move(genome);
        genomeData[targetIndex] = DynamicGenomeData{};  // Placeholder
        
        // Create actual metadata with placed genome
        DynamicGenomeData metadata = metadataCreator(genomes[targetIndex], targetIndex);
        
        // Update with actual metadata
        genomeData[targetIndex] = std::move(metadata);
        
        // Perform fitness evaluation if genome is evaluable
        if (!genomeData[targetIndex].isUnderRepair) {
            FitnessResultType fitness = fitnessStrategy->evaluate(
                genomes[targetIndex].get_phenotype(), 
                speciationControl
            );
            
            // Insert fitness result into multimap
            auto& fitnessResults = container.getFitnessResults(currentGeneration);
            fitnessResults.insert({fitness, targetIndex});
        }
    }
    
    LOG_DEBUG("genomePlacement: Successfully placed genome at index {}", targetIndex);
}

}