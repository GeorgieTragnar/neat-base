#pragma once

#include <functional>
#include <cstdint>
#include <cassert>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/PopulationData.hpp"

namespace Operator {

template<typename FitnessResultType>
void mutationPlacement(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& currentGeneration,
    const std::function<Genome(const Genome&)>& mutationFunction,
    const Genome& parentGenome,
    const DynamicGenomeData& parentData
);

template<typename FitnessResultType>
void mutationPlacement(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& currentGeneration,
    const std::function<Genome(const Genome&)>& mutationFunction,
    const Genome& parentGenome,
    const DynamicGenomeData& parentData
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
    
    // Execute the mutation function on parent genome
    Genome offspring = mutationFunction(parentGenome);
    
    // Create offspring metadata based on parent data
    DynamicGenomeData offspringData = parentData;
    offspringData.parentAIndex = static_cast<uint32_t>(targetIndex);  // Same index in last generation
    offspringData.parentBIndex = UINT32_MAX;                          // Single parent (not crossover)
    // Note: isUnderRepair and other flags will be updated by other operators
    
    // Place offspring and metadata directly into container at target index
    currentGenomes[targetIndex] = std::move(offspring);
    currentGenomeData[targetIndex] = offspringData;
}

} // namespace Operator