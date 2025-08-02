#pragma once

#include <functional>
#include <cstdint>
#include <cassert>
#include <memory>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/PopulationData.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/analysis/FitnessStrategy.hpp"
#include "version3/analysis/SpeciationControlUnit.hpp"

namespace Operator {

template<typename FitnessResultType>
void mutationPlacement(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& currentGeneration,
    Genome offspring,
    const DynamicGenomeData& parentData,
    std::function<void(Genome&, DynamicGenomeData&)> postPlacementOperator,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl,
    GlobalIndexRegistry& registry
);

template<typename FitnessResultType>
void mutationPlacement(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& currentGeneration,
    Genome offspring,
    const DynamicGenomeData& parentData,
    std::function<void(Genome&, DynamicGenomeData&)> postPlacementOperator,
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    Analysis::SpeciationControlUnit& speciationControl,
    GlobalIndexRegistry& registry
) {
    // Validate generation bounds
    assert(currentGeneration > 0 && 
           "mutationPlacement: Cannot place in generation 0 during evolution");
    
    // Get mutable access to current generation
    auto& currentGenomes = container.getCurrentGenomes(currentGeneration);
    auto& currentGenomeData = container.getCurrentGenomeData(currentGeneration);
    
    // Validate bounds
    assert(targetIndex < currentGenomes.size() && 
           "mutationPlacement: targetIndex out of bounds");
    assert(currentGenomes.size() == currentGenomeData.size() && 
           "mutationPlacement: Genome and data vectors must have equal size");
    
    // Validate parent data has proper genomeIndex before copying
    assert(parentData.genomeIndex != UINT32_MAX && "Parent genome must have valid genomeIndex before mutation placement");
    
    // Copy parent data and place genome with move semantics
    currentGenomeData[targetIndex] = parentData;
    currentGenomes[targetIndex] = std::move(offspring);
    
    // MutationPlacement operator handles parent index setup
    currentGenomeData[targetIndex].parentAIndex = static_cast<uint32_t>(targetIndex);  // Same index in last generation
    currentGenomeData[targetIndex].parentBIndex = UINT32_MAX;                          // Single parent (not crossover)
    
    // Store parent species before post-placement operations
    uint32_t parentSpeciesId = parentData.speciesId;
    
    // Execute post-placement operations (species assignment, cycle detection, phenotype updates, etc.)
    postPlacementOperator(
        currentGenomes[targetIndex],
        currentGenomeData[targetIndex]
    );
    
    // Reset elimination counter if offspring moved to a different species
    if (currentGenomeData[targetIndex].speciesId != parentSpeciesId) {
        currentGenomeData[targetIndex].pendingEliminationCounter = 0;
    }
    
    // Perform fitness evaluation if genome is evaluable
    if (!currentGenomeData[targetIndex].isUnderRepair) {
        if (currentGenomeData[targetIndex].isElite) {
            // Elite genomes: use previous generation's fitness
            const auto& lastFitnessResults = container.getLastFitnessResults(currentGeneration);
            FitnessResultType eliteFitness{};
            bool fitnessFound = false;
            
            for (const auto& [fitness, globalIndex] : lastFitnessResults) {
                if (globalIndex == targetIndex) {
                    eliteFitness = fitness;
                    fitnessFound = true;
                    break;
                }
            }
            assert(fitnessFound && "Elite genome must have fitness from previous generation");
            
            auto& fitnessResults = container.getCurrentFitnessResults(currentGeneration);
            fitnessResults.insert({eliteFitness, targetIndex});
        } else {
            // Non-elite genomes: evaluate fitness
            FitnessResultType fitness = fitnessStrategy->evaluate(
                currentGenomes[targetIndex].get_phenotype(), 
                speciationControl
            );
            
            auto& fitnessResults = container.getCurrentFitnessResults(currentGeneration);
            fitnessResults.insert({fitness, targetIndex});
        }
    }
}

} // namespace Operator