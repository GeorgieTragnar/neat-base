#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>
#include "../data/PopulationData.hpp"
#include "../data/PopulationContainer.hpp"
#include "../data/GlobalIndexRegistry.hpp"
#include "GenerationPopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

// Parameters for species equilibrium elimination
struct SpeciesEquilibriumParams {
    uint32_t equilibriumSpeciesCount;                    // Target number of species
    double speciesPendingEliminationPercentage;          // Percentage of excess species to target
    double elitePlacementProtectionPercentage;           // Champion protection threshold
    uint32_t maxPendingEliminationRating;               // Max rating before elimination
    uint32_t targetPopulationSize;                      // Target population size
    
    SpeciesEquilibriumParams() = delete;
    
    SpeciesEquilibriumParams(
        uint32_t equilibriumSpeciesCount,
        double speciesPendingEliminationPercentage,
        double elitePlacementProtectionPercentage,
        uint32_t maxPendingEliminationRating,
        uint32_t targetPopulationSize
    ) : equilibriumSpeciesCount(equilibriumSpeciesCount),
        speciesPendingEliminationPercentage(speciesPendingEliminationPercentage),
        elitePlacementProtectionPercentage(elitePlacementProtectionPercentage),
        maxPendingEliminationRating(maxPendingEliminationRating),
        targetPopulationSize(targetPopulationSize) {
        
        assert(elitePlacementProtectionPercentage >= 0.0 && elitePlacementProtectionPercentage <= 1.0);
        assert(speciesPendingEliminationPercentage >= 0.0 && speciesPendingEliminationPercentage <= 1.0);
        assert(equilibriumSpeciesCount > 0);
        assert(targetPopulationSize > 0);
    }
};

// SpeciesEquilibriumElimination operator: Eliminates worst-performing species based on equilibrium
template<typename FitnessResultType>
void speciesEquilibriumElimination(
    const GenerationPopulationData<FitnessResultType>& populationData,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const SpeciesEquilibriumParams& params,
    GlobalIndexRegistry& registry
);

// Template implementation
template<typename FitnessResultType>
void speciesEquilibriumElimination(
    const GenerationPopulationData<FitnessResultType>& populationData,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    const SpeciesEquilibriumParams& params,
    GlobalIndexRegistry& registry
) {
    auto logger = LOGGER("population.SpeciesEquilibriumElimination");
    
    // Get genome data from container
    const auto& genomeData = container.getGenomeData(generation);
    
    // Calculate excess species count vs equilibrium target
    const uint32_t excessSpeciesCount = (populationData.activeSpeciesCount > params.equilibriumSpeciesCount)
        ? populationData.activeSpeciesCount - params.equilibriumSpeciesCount : 0;
    
    assert(excessSpeciesCount <= populationData.activeSpeciesCount && 
           "there cannot be more excess than active species");
    
    const uint32_t speciesPendingElimination = excessSpeciesCount > 0
        ? static_cast<uint32_t>(excessSpeciesCount * params.speciesPendingEliminationPercentage) : 0;
    
    assert(speciesPendingElimination <= populationData.activeSpeciesCount && 
           "there cannot be more species pending elimination than active species");
    
    // Track species that had their rating increased this generation
    std::unordered_set<uint32_t> speciesWithIncreasedRating;
    std::unordered_set<uint32_t> eliminatedSpeciesIds;
    
    if (speciesPendingElimination > 0) {
        // Create vector of active species sorted by average rank (worst to best performance)
        std::vector<std::pair<uint32_t, double>> speciesRanking;
        for (const auto& [speciesId, avgRank] : populationData.speciesAverageRanks) {
            // Only consider active species that aren't already marked for elimination
            auto speciesIt = speciesData.find(speciesId);
            if (speciesIt != speciesData.end() && !speciesIt->second.isMarkedForElimination && 
                speciesIt->second.currentPopulationSize > 0) {
                speciesRanking.emplace_back(speciesId, avgRank);
            }
        }
        
        // Sort by average rank (ascending - worst performers first)
        std::sort(speciesRanking.begin(), speciesRanking.end(), 
                  [](const auto& a, const auto& b) { return a.second < b.second; });
        
        // Check the worst performing species for elite placement protection
        const size_t totalPopulation = populationData.fitnessResults.size();
        const size_t speciesToCheck = std::min(static_cast<size_t>(speciesPendingElimination), 
                                               speciesRanking.size());
        
        for (size_t i = 0; i < speciesToCheck; ++i) {
            const uint32_t speciesId = speciesRanking[i].first;
            
            // Get the best genome rank for this species (their champion)
            auto championRankIt = populationData.speciesBestGenomeRank.find(speciesId);
            if (championRankIt == populationData.speciesBestGenomeRank.end()) {
                // This shouldn't happen, but handle gracefully
                LOG_WARN("Species {} has no champion rank recorded, skipping elimination check", speciesId);
                continue;
            }
            
            const size_t championRank = championRankIt->second;
            
            // Calculate champion's percentile in population (0.0 = worst, 1.0 = best)
            const double championPercentile = (totalPopulation > 1) ? 
                static_cast<double>(championRank) / (totalPopulation - 1) : 1.0;
            
            // Check if champion is protected by elite placement percentage
            if (championPercentile >= params.elitePlacementProtectionPercentage) {
                // Champion is in top X% of population - species is protected from elimination
                LOG_DEBUG("Species {} protected: champion at rank {} (percentile {:.1f}%) >= protection threshold {:.1f}%", 
                         speciesId, championRank, championPercentile * 100.0, 
                         params.elitePlacementProtectionPercentage * 100.0);
                continue;
            }
            
            // Species is not protected - increase pending elimination rating
            auto& speciesInfo = speciesData[speciesId];
            speciesInfo.pendingEliminationRating++;
            speciesWithIncreasedRating.insert(speciesId);
            
            LOG_DEBUG("Species {} pending elimination increased to {}: champion at rank {} (percentile {:.1f}%) < protection threshold {:.1f}%", 
                     speciesId, speciesInfo.pendingEliminationRating, championRank, 
                     championPercentile * 100.0, params.elitePlacementProtectionPercentage * 100.0);
            
            // Mark species for elimination if rating exceeds limit
            if (speciesInfo.pendingEliminationRating >= params.maxPendingEliminationRating) {
                // Check if species being eliminated contains Elite genomes
                bool hasEliteGenomes = false;
                for (const auto& [fitnessResult, globalIndex] : populationData.fitnessResults) {
                    if (genomeData[globalIndex].speciesId == speciesId && 
                        genomeData[globalIndex].isElite) {
                        hasEliteGenomes = true;
                        LOG_DEBUG("ELITE_TRACK: SpeciesEquilibriumElimination - Species {} has Elite genome {} but being eliminated!", 
                                 speciesId, globalIndex);
                    }
                }
                
                speciesInfo.isMarkedForElimination = true;
                eliminatedSpeciesIds.insert(speciesId);
                LOG_INFO("Species {} marked for elimination: pending rating {} >= limit {} (hasElites: {})", 
                        speciesId, speciesInfo.pendingEliminationRating, params.maxPendingEliminationRating, hasEliteGenomes);
            }
        }
    }
    
    // Decrease pending elimination rating for active species that weren't targeted for elimination
    for (auto& [speciesId, data] : speciesData) {
        // Only process active species that aren't marked for elimination and didn't have rating increased
        if (data.currentPopulationSize > 0 && !data.isMarkedForElimination && 
            speciesWithIncreasedRating.find(speciesId) == speciesWithIncreasedRating.end()) {
            
            if (data.pendingEliminationRating > 0) {
                data.pendingEliminationRating--;
                LOG_DEBUG("Species {} pending elimination decreased to {}: not targeted for elimination this generation", 
                         speciesId, data.pendingEliminationRating);
            }
        }
    }
}

} // namespace Operator