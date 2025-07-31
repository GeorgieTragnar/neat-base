#pragma once

#include <functional>
#include <memory>
#include <cassert>
#include "../data/Genome.hpp"
#include "../data/PopulationContainer.hpp"
#include "../data/GlobalIndexRegistry.hpp"
#include "../data/HistoryTracker.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

template<typename FitnessResultType>
void bootstrap(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    uint32_t& generation,
    uint32_t targetPopulationSize,
    std::function<Genome()> createGenome,
    std::function<DynamicGenomeData(const Genome&, uint32_t)> createMetadata
);

template<typename FitnessResultType>
inline void bootstrap(
    PopulationContainer<FitnessResultType>& container,
    GlobalIndexRegistry& registry,
    uint32_t& generation,
    uint32_t targetPopulationSize,
    std::function<Genome()> createGenome,
    std::function<DynamicGenomeData(const Genome&, uint32_t)> createMetadata
) {
    auto logger = LOGGER("population.Bootstrap");
    
    assert(targetPopulationSize > 0 && "Bootstrap requires positive population size");
    assert(createGenome != nullptr && "Bootstrap requires valid genome creator function");
    assert(createMetadata != nullptr && "Bootstrap requires valid metadata creator function");
    
    LOG_INFO("Bootstrap: Creating initial population of {} genomes", targetPopulationSize);
    
    // Generation 0: Create initial population
    // Constructor already initialized all arrays to zero - perfect for gen 0
    LOG_DEBUG("Bootstrap: Creating generation 0");
    container.clearGenerationFitnessResults(0);
    
    for (uint32_t i = 0; i < targetPopulationSize; ++i) {
        // Create genome via user-provided lambda
        Genome genome = createGenome();
        
        // Add to container with placeholder metadata (returns global index)
        DynamicGenomeData placeholderMetadata{};
        uint32_t globalIndex = container.push_back(0, std::move(genome), std::move(placeholderMetadata));
        
        // Verify index matches expectation (sequential for bootstrap)
        assert(globalIndex == i && "Bootstrap expects sequential global indices");
        
        // Get reference to placed genome
        const auto& genomes = container.getGenomes(0);
        
        // Create actual metadata via user-provided lambda
        DynamicGenomeData metadata = createMetadata(genomes[globalIndex], globalIndex);
        
        // Update with actual metadata
        auto& genomeData = container.getGenomeData(0);
        genomeData[globalIndex] = std::move(metadata);
        
        // Ensure genome is marked as Active in registry
        assert(registry.getState(globalIndex) == GenomeState::Active && 
               "Bootstrap genomes must be Active");
    }
    
    LOG_DEBUG("Bootstrap: Generation 0 complete with {} genomes", targetPopulationSize);
    
    // PopulationContainer automatically synchronized all 3 generations during push_back operations
    // Generation 1 is ready for announceNewGeneration(1) call
    generation = 1;
    LOG_INFO("Bootstrap: Complete. Evolution will start at generation {}", generation);
    
    // Validate container state - all generations should have same population due to auto-sync
    assert(container.getGenerationSize(0) == targetPopulationSize);
    assert(container.getGenerationSize(1) == targetPopulationSize);
    assert(container.getGenerationSize(2) == targetPopulationSize);
    assert(container.getFitnessResults(0).empty() && "Generation 0 should have no fitness results");
    assert(container.getFitnessResults(1).empty() && "Generation 1 should have no fitness results");
    assert(container.getFitnessResults(2).empty() && "Generation 2 should have no fitness results");
}

} // namespace Operator