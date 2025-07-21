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
    Genome offspring,
    const DynamicGenomeData& parentData,
    std::function<void(Genome&, DynamicGenomeData&)> postPlacementOperator
);

template<typename FitnessResultType>
void mutationPlacement(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& currentGeneration,
    Genome offspring,
    const DynamicGenomeData& parentData,
    std::function<void(Genome&, DynamicGenomeData&)> postPlacementOperator
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
    
    // Copy parent data and place genome with move semantics
    currentGenomeData[targetIndex] = parentData;
    currentGenomes[targetIndex] = std::move(offspring);
    
    // MutationPlacement operator handles parent index setup
    currentGenomeData[targetIndex].parentAIndex = static_cast<uint32_t>(targetIndex);  // Same index in last generation
    currentGenomeData[targetIndex].parentBIndex = UINT32_MAX;                          // Single parent (not crossover)
    
    // Execute post-placement operations (species assignment, cycle detection, phenotype updates, etc.)
    postPlacementOperator(
        currentGenomes[targetIndex],
        currentGenomeData[targetIndex]
    );
}

} // namespace Operator