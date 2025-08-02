#pragma once

#include <unordered_map>
#include <vector>
#include <cassert>
#include "../data/PopulationData.hpp"
#include "../data/PopulationContainer.hpp"
#include "../data/GlobalIndexRegistry.hpp"
#include "GenerationPopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

// Parameters for genome equilibrium elimination
struct GenomeEquilibriumParams {
    uint32_t targetPopulationSize;                      // Target population size
    double genomesPendingEliminationPercentage;         // Percentage of excess genomes to target
    uint32_t maxPendingEliminationCounter;             // Max counter before elimination
    
    GenomeEquilibriumParams() = delete;
    
    GenomeEquilibriumParams(
        uint32_t targetPopulationSize,
        double genomesPendingEliminationPercentage,
        uint32_t maxPendingEliminationCounter
    ) : targetPopulationSize(targetPopulationSize),
        genomesPendingEliminationPercentage(genomesPendingEliminationPercentage),
        maxPendingEliminationCounter(maxPendingEliminationCounter) {
        
        assert(genomesPendingEliminationPercentage >= 0.0 && genomesPendingEliminationPercentage <= 1.0);
        assert(targetPopulationSize > 0);
    }
};

// GenomeEquilibriumElimination operator: Eliminates worst-performing genomes within species
template<typename FitnessResultType>
void genomeEquilibriumElimination(
    const GenerationPopulationData<FitnessResultType>& populationData,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const GenomeEquilibriumParams& params,
    GlobalIndexRegistry& registry
);

// Template implementation
template<typename FitnessResultType>
void genomeEquilibriumElimination(
    const GenerationPopulationData<FitnessResultType>& populationData,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const GenomeEquilibriumParams& params,
    GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("population.GenomeEquilibriumElimination");
    
    // Get genome data from container (need mutable access for counter updates)
    auto& genomeData = container.getCurrentGenomeData(generation);
    
    // Calculate per-species equilibrium (target population / active species count)
    const uint32_t equilibrium = (populationData.activeSpeciesCount == 0) ? 0 :
        (populationData.activeSpeciesCount > params.targetPopulationSize) ? 1 : 
        params.targetPopulationSize / populationData.activeSpeciesCount;
    
    // Process each active species for genome-level elimination
    for (const auto& [speciesId, speciesInfo] : speciesData) {
        // Skip extinct, marked for elimination, or species without valid genomes
        if (speciesInfo.currentPopulationSize == 0 || 
            populationData.speciesGrouping.find(speciesId) == populationData.speciesGrouping.end()) {
            continue;
        }
        
        if (speciesInfo.isMarkedForElimination) {
            auto itSpeciesInfo = populationData.speciesGrouping.find(speciesId);
            if (itSpeciesInfo != populationData.speciesGrouping.end()) {
                for (auto& genomeIndex : itSpeciesInfo->second) {
                    registry.markForElimination(genomeIndex);
                }
            }
            continue;
        }

        const std::vector<uint32_t>& validGenomes = populationData.speciesGrouping.at(speciesId);
        const uint32_t activeGenomeCount = static_cast<uint32_t>(validGenomes.size());
        const uint32_t excessGenomeCount = (activeGenomeCount > equilibrium) 
            ? activeGenomeCount - equilibrium : 0;
        
        if (excessGenomeCount == 0) {
            // No excess genomes - still need to decrease pending counters for all genomes
            for (size_t i = 0; i < validGenomes.size(); ++i) {
                if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                    continue;
                }
                
                auto& genome = genomeData[validGenomes[i]];
                assert(!genome.isUnderRepair && "invalid access to genome data under repair");
                
                if (genome.pendingEliminationCounter > 0) {
                    genome.pendingEliminationCounter--;
                }
            }
            continue;
        }
        
        assert(excessGenomeCount < activeGenomeCount && 
               "we cannot have more excess genomes than the active species population size");
        
        const uint32_t genomesPendingElimination = static_cast<uint32_t>(
            excessGenomeCount * params.genomesPendingEliminationPercentage);
        
        if (genomesPendingElimination == 0) {
            // No genomes targeted for pending elimination - decrease all counters
            for (size_t i = 0; i < validGenomes.size(); ++i) {
                if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                    continue;
                }
                
                auto& genome = genomeData[validGenomes[i]];
                assert(!genome.isUnderRepair && "invalid access to genome data under repair");
                
                if (genome.pendingEliminationCounter > 0) {
                    genome.pendingEliminationCounter--;
                }
            }
            continue;
        }
        
        assert(genomesPendingElimination < activeGenomeCount && 
               "we cannot mark more genomes for pending elimination than the active species population size");
        
        assert(!validGenomes.empty() && 
               "at this point species without valid genomes should have been skipped already");
        
        LOG_DEBUG("Species {}: {} genomes, equilibrium {}, excess {}, targeting {} for pending elimination",
                 speciesId, activeGenomeCount, equilibrium, excessGenomeCount, genomesPendingElimination);
        
        // Increment pending elimination counter for worst performers (beginning of fitness-ordered vector)
        for (size_t i = 0; i < genomesPendingElimination; ++i) {
            // Skip genomes that may have been marked for elimination during this update cycle
            if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                continue;
            }
            
            auto& genome = genomeData[validGenomes[i]];
            assert(!genome.isUnderRepair && "invalid access to genome data under repair");
            
            genome.pendingEliminationCounter++;
            
            // Mark for elimination if counter exceeds limit
            if (genome.pendingEliminationCounter >= params.maxPendingEliminationCounter) {
                registry.markForElimination(validGenomes[i]);
                LOG_DEBUG("Genome {} in species {} marked for elimination: counter {} >= limit {}",
                         validGenomes[i], speciesId, genome.pendingEliminationCounter, 
                         params.maxPendingEliminationCounter);
            }
        }
        
        // Decrement pending elimination counter for better performers (end of fitness-ordered vector)
        for (size_t i = genomesPendingElimination; i < validGenomes.size(); ++i) {
            // Skip genomes that may have been marked for elimination during this update cycle
            if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                continue;
            }
            
            auto& genome = genomeData[validGenomes[i]];
            assert(!genome.isUnderRepair && "invalid access to genome data under repair");
            
            if (genome.pendingEliminationCounter > 0) {
                genome.pendingEliminationCounter--;
            }
        }
    }
}

} // namespace Operator