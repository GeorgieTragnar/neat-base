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
    std::function<DynamicGenomeData(const Genome&, size_t)> metadataCreator,
    uint32_t currentGeneration,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl
);

template<typename FitnessResultType>
inline void genomePlacement(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    Genome genome,
    std::function<DynamicGenomeData(const Genome&, size_t)> metadataCreator,
    uint32_t currentGeneration,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl
) {
    auto logger = LOGGER("operator.GenomePlacement");
    
    // First, try to get a free index for replacement
    uint32_t targetIndex = registry.getFreeIndex();
    size_t actualIndex;
    
    if (targetIndex == INVALID_INDEX) {
        // No free indices available - add new genome to population
        actualIndex = container.getGenomes(currentGeneration).size();
        LOG_DEBUG("genomePlacement: Adding new genome to population at index {}", actualIndex);
        
        // Add genome to container first with placeholder metadata
        DynamicGenomeData placeholderMetadata{};
        container.push_back(currentGeneration, std::move(genome), std::move(placeholderMetadata));
        
        // Get reference to placed genome and create actual metadata
        auto& genomes = container.getGenomes(currentGeneration);
        auto& genomeData = container.getGenomeData(currentGeneration);
        DynamicGenomeData metadata = metadataCreator(genomes[actualIndex], actualIndex);
        
        // Update with actual metadata
        genomeData[actualIndex] = std::move(metadata);
        
        // Perform fitness evaluation if genome is evaluable
        if (!genomeData[actualIndex].isUnderRepair) {
            FitnessResultType fitness = fitnessStrategy->evaluate(
                genomes[actualIndex].get_phenotype(), 
                speciationControl
            );
            
            // Insert fitness result into multimap
            auto& fitnessResults = container.getFitnessResults(currentGeneration);
            fitnessResults.insert({fitness, actualIndex});
        }
        
    } else {
        // Using existing ReadyForReplacement slot
        actualIndex = targetIndex;
        LOG_DEBUG("genomePlacement: Replacing genome at index {} (reusing ReadyForReplacement slot)", actualIndex);
        
        // Verify the index is actually ready for replacement
        GenomeState state = registry.getState(targetIndex);
        assert(state == GenomeState::ReadyForReplacement && 
               "genomePlacement: Attempted to replace genome that is not marked ReadyForReplacement");
        
        // Replace genome in container first with placeholder metadata
        auto& genomes = container.getGenomes(currentGeneration);
        auto& genomeData = container.getGenomeData(currentGeneration);
        genomes[targetIndex] = std::move(genome);
        genomeData[targetIndex] = DynamicGenomeData{};  // Placeholder
        
        // Create actual metadata with placed genome
        DynamicGenomeData metadata = metadataCreator(genomes[targetIndex], actualIndex);
        
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
    
    LOG_DEBUG("genomePlacement: Successfully placed genome at index {}", actualIndex);
}

}