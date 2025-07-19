#pragma once

#include <cstdint>
#include <limits>

#include "version3/data/Genome.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/PopulationData.hpp"
#include "version3/evolution/Crossover.hpp"
#include "version3/validation/CycleDetection.hpp"
#include "version3/phenotype/PhenotypeConstruct.hpp"

namespace Operator {

// Forward declarations
class CrossoverManagementParams;

template<typename FitnessResultType>
void crossoverManagement(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& currentGeneration,
    const CrossoverManagementParams& params
);

// Parameters class definition
class CrossoverManagementParams {
public:
    CrossoverManagementParams(
        const CrossoverParams& crossoverParams,
        const CycleDetectionParams& cycleParams
    ) : _crossoverParams(crossoverParams),
        _cycleParams(cycleParams) {}

    CrossoverManagementParams(const CrossoverManagementParams& other) = default;
    CrossoverManagementParams& operator=(const CrossoverManagementParams& other) = default;

protected:
    const CrossoverParams _crossoverParams;
    const CycleDetectionParams _cycleParams;

    template<typename FitnessResultType>
    friend void crossoverManagement(
        PopulationContainer<FitnessResultType>& container,
        GlobalIndexRegistry& registry,
        const size_t& parentAIndex,
        const size_t& parentBIndex,
        const uint32_t& currentGeneration,
        const CrossoverManagementParams& params
    );
};

// Template implementation
template<typename FitnessResultType>
void crossoverManagement(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& currentGeneration,
    const CrossoverManagementParams& params
) {
    // Access current generation data
    const auto& genomes = container.getGenomes(currentGeneration);
    const auto& genomeData = container.getGenomeData(currentGeneration);
    const auto& fitnessResults = container.getFitnessResults(currentGeneration);
    
    // Validate parent indices
    assert(parentAIndex < genomes.size() && "Parent A index out of bounds");
    assert(parentBIndex < genomes.size() && "Parent B index out of bounds");
    
    // Validate parent states
    auto parentAState = registry.getState(parentAIndex);
    auto parentBState = registry.getState(parentBIndex);
    assert((parentAState == GenomeState::Active || parentAState == GenomeState::Elite) && 
           "Crossover parent A should be in Active or Elite state");
    assert((parentBState == GenomeState::Active || parentBState == GenomeState::Elite) && 
           "Crossover parent B should be in Active or Elite state");
    
    // Get parent genomes and data
    const Genome& parentA = genomes[parentAIndex];
    const Genome& parentB = genomes[parentBIndex];
    const DynamicGenomeData& parentAData = genomeData[parentAIndex];
    const DynamicGenomeData& parentBData = genomeData[parentBIndex];
    
    // Extract fitness values from multimap
    FitnessResultType fitnessA{}, fitnessB{};
    bool foundA = false, foundB = false;
    for (const auto& [fitness, globalIndex] : fitnessResults) {
        if (globalIndex == parentAIndex) {
            fitnessA = fitness;
            foundA = true;
        }
        if (globalIndex == parentBIndex) {
            fitnessB = fitness;
            foundB = true;
        }
        if (foundA && foundB) break;
    }
    assert(foundA && "Parent A fitness not found in results");
    assert(foundB && "Parent B fitness not found in results");
    
    // Perform crossover
    Genome offspring = Operator::crossover(
        parentA, fitnessA,
        parentB, fitnessB,
        params._crossoverParams
    );
    
    // Check for cycles
    bool hasCycles = Operator::hasCycles(offspring, params._cycleParams);
    
    // Construct phenotype if no cycles
    if (!hasCycles) {
        Operator::phenotypeConstruct(offspring);
    }
    
    // Create offspring metadata
    DynamicGenomeData offspringData;
    offspringData.speciesId = parentAData.speciesId;  // Inherit from parent A
    offspringData.pendingEliminationCounter = 0;      // Fresh start
    offspringData.isUnderRepair = hasCycles;
    offspringData.isMarkedForElimination = false;
    offspringData.parentAIndex = parentAIndex;
    offspringData.parentBIndex = parentBIndex;
    
    // Get target index for placement
    uint32_t targetIndex = registry.getFreeIndex();
    
    if (targetIndex == INVALID_INDEX) {
        // No free indices - add new genome to population
        container.push_back(currentGeneration, std::move(offspring), std::move(offspringData));
    } else {
        // Reuse existing slot (already marked Active by getFreeIndex)
        auto& mutableGenomes = container.getGenomes(currentGeneration);
        auto& mutableGenomeData = container.getGenomeData(currentGeneration);
        mutableGenomes[targetIndex] = std::move(offspring);
        mutableGenomeData[targetIndex] = offspringData;
    }
}

} // namespace Operator