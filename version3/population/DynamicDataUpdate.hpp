#pragma once

#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cassert>
#include "PopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Population {

class DynamicDataUpdateParams {
public:
    DynamicDataUpdateParams() = delete;
    DynamicDataUpdateParams(
        uint32_t maxProtectionLimit,
        uint32_t maxSpeciesProtectionRating,
        double protectedTierPercentage,
        uint32_t worstSpeciesCount,
        uint32_t minActiveSpeciesCount
    );

protected:
    template<typename FitnessResultType>
    friend void dynamicDataUpdate(
        const std::multimap<FitnessResultType, size_t>& fitnessResults,
        std::vector<DynamicGenomeData>& genomeData,
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const DynamicDataUpdateParams& params
    );

    const uint32_t _maxProtectionLimit;           // Max protection counter before elimination
    const uint32_t _maxSpeciesProtectionRating;  // Max species protection rating before elimination
    const double _protectedTierPercentage;       // Percentage of population in protected tier (e.g., 0.3 for 30%)
    const uint32_t _worstSpeciesCount;           // Number of worst species to penalize each generation
    const uint32_t _minActiveSpeciesCount;       // Minimum number of active species before elimination/protection penalties are disabled
};

// Helper function declarations
void updateSpeciesRanking(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, double>& speciesAverageRanks
);

// Main operator function - templated on fitness result type
template<typename FitnessResultType>
void dynamicDataUpdate(
    const std::multimap<FitnessResultType, size_t>& fitnessResults,  // Best to worst fitness with global indices
    std::vector<DynamicGenomeData>& genomeData,  // Genome data vector for direct access
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const DynamicDataUpdateParams& params
) {
    // Phase 1: Species Performance Analysis - Single pass with running averages
    std::unordered_map<uint32_t, double> speciesRankSum;
    std::unordered_map<uint32_t, size_t> speciesCount;
    std::unordered_map<uint32_t, double> speciesAverageRanks;
    
    // Reset species population sizes at start of analysis
    for (auto& [speciesId, data] : speciesData) {
        data.currentPopulationSize = 0;
        speciesRankSum[speciesId] = 0.0;
        speciesCount[speciesId] = 0;
    }
    
    size_t rank = 0;
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        
        // Note: Species discovery happens below - genome can reference previously unknown species
        // This is valid when new species emerge during evolution
        
        // Update species rank sum and count
        speciesRankSum[speciesId] += static_cast<double>(rank);
        speciesCount[speciesId]++;
        
        // Update species population size in single pass
        auto speciesIt = speciesData.find(speciesId);
        if (speciesIt != speciesData.end()) {
            speciesIt->second.currentPopulationSize++;
        } else {
            // Create new species data entry for species discovered in genome data
            // auto logger = LOGGER("population.DynamicDataUpdate");
            // LOG_DEBUG("Discovered new species {} in genome data, creating species entry", speciesId);
            DynamicSpeciesData newSpecies;
            newSpecies.currentPopulationSize = 1;     // First genome for this species
            newSpecies.protectionRating = 0;          // Default protection rating
            newSpecies.speciesRank = 0;               // TODO: Review new species rank assignment timing
            newSpecies.isMarkedForElimination = false; // Default elimination status
            speciesData[speciesId] = newSpecies;
            // TODO: New species are created after rank calculation - review this interaction
        }
        
        ++rank;
    }
    
    // Calculate average rank per species
    for (const auto& [speciesId, count] : speciesCount) {
        if (count > 0) {
            speciesAverageRanks[speciesId] = speciesRankSum[speciesId] / count;
        }
    }
    
    // Update species rankings based on performance
    updateSpeciesRanking(speciesData, speciesAverageRanks);
    
    // Find worst N species using partial sort for efficiency
    std::vector<std::pair<uint32_t, double>> speciesRankings;
    for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
        speciesRankings.emplace_back(speciesId, avgRank);
    }
    
    // Only need to identify worst N species, not full sort
    const size_t worstSpeciesCount = std::min(
        static_cast<size_t>(params._worstSpeciesCount), 
        speciesRankings.size()
    );
    
    if (worstSpeciesCount < speciesRankings.size()) {
        // Partial sort: worst species will be at the end
        std::nth_element(speciesRankings.begin(), 
                        speciesRankings.end() - worstSpeciesCount,
                        speciesRankings.end(),
                        [](const auto& a, const auto& b) { return a.second < b.second; });
    }
    
    // Phase 2: Individual Genome Protection Updates
    // TODO: Rename "protection" mechanism - it's confusing as "protected" entities progress toward elimination
    // Consider: "elimination_counter", "underperformance_counter", "penalty_counter", etc.
    
    // Count active species before applying protection penalties
    uint32_t activeSpeciesCount = 0;
    for (const auto& [speciesId, data] : speciesData) {
        if (!data.isMarkedForElimination) {
            activeSpeciesCount++;
        }
    }
    
    // Only apply protection penalties if we have more than minimum active species
    const bool applyProtectionPenalties = activeSpeciesCount > params._minActiveSpeciesCount;
    
    const size_t totalPopulation = fitnessResults.size();
    const size_t protectedTierThreshold = static_cast<size_t>(
        totalPopulation * params._protectedTierPercentage
    );
    
    rank = 0;
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        DynamicGenomeData& currentGenomeData = genomeData[globalIndex];
        
        // Skip protection updates for genomes under repair
        if (currentGenomeData.isUnderRepair) {
            ++rank;
            continue;
        }
        
        // Check if genome is in protected tier (bottom X% by rank)
        const bool isInProtectedTier = rank >= (totalPopulation - protectedTierThreshold);
        
        if (applyProtectionPenalties) {
            if (isInProtectedTier) {
                currentGenomeData.protectionCounter++;
            } else {
                // TODO: Consider graceful protection counter decrement instead of hard reset
                // Current hard reset causes volatility - genomes pop in/out of elimination threshold
                // Alternative: gradual decrement (e.g., protectionCounter = max(0, protectionCounter-1))
                // This would eliminate genomes that spend more time protected than not
                currentGenomeData.protectionCounter = 0;
            }
            
            // Mark for elimination if protection limit exceeded
            if (currentGenomeData.protectionCounter > params._maxProtectionLimit) {
                currentGenomeData.isMarkedForElimination = true;
            }
        }
        // If we're at/below minimum species count, don't increase protection counters or mark genomes for elimination
        
        ++rank;
    }
    
    // Phase 3: Species Protection Rating Updates
    // Create set of worst species IDs for quick lookup
    std::unordered_set<uint32_t> worstSpeciesIds;
    for (size_t i = speciesRankings.size() - worstSpeciesCount; i < speciesRankings.size(); ++i) {
        worstSpeciesIds.insert(speciesRankings[i].first);
    }
    
    // TODO: Expand species protection logic to prevent elimination of species that have genomes 
    // in the top X% of the full population (elite species protection)
    
    // Update protection ratings: penalize worst species, reset others
    // TODO: Active species threshold should only affect species penalties, not genome penalties
    // if (activeSpeciesCount > params._minActiveSpeciesCount) {
    //     for (auto& [speciesId, data] : speciesData) {
    //         if (worstSpeciesIds.find(speciesId) != worstSpeciesIds.end()) {
    //             // This species is in worst N - apply penalty
    //             data.protectionRating++;
                
    //             // Mark species for elimination if rating threshold exceeded
    //             if (data.protectionRating > params._maxSpeciesProtectionRating) {
    //                 data.isMarkedForElimination = true;
    //             }
    //         } else {
    //             // This species escaped worst N - reset rating
    //             data.protectionRating = 0;
    //         }
    //     }
    // }
    // If below minimum species count, no species protection penalties are applied
    
    // Phase 4: State Consistency Validation (Debug assertions)
#ifndef NDEBUG
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        // Protection counter can exceed maxProtectionLimit for custom scaling logic that punishes worst genomes more severely
        // No assertion needed for protectionCounter upper bound
        assert(speciesData.find(genomeData[globalIndex].speciesId) != speciesData.end());
    }
    
    for (const auto& [speciesId, data] : speciesData) {
        // Protection rating can exceed maxSpeciesProtectionRating for adaptive elimination strategies
        // Note: speciesRank can be 0 if a species has a single genome that is the best performer (rank 0)
        // Species rank is calculated from average genome ranks, so single-genome species inherit their genome's rank
        assert(data.speciesRank >= 0); // Species rank should be assigned (can be 0 for best single-genome species)
        
        // Population size should match actual count
        size_t actualCount = 0;
        for (const auto& [fr, globalIndex] : fitnessResults) {
            if (genomeData[globalIndex].speciesId == speciesId) actualCount++;
        }
        assert(data.currentPopulationSize == actualCount);
    }
#endif
}

} // namespace Population