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

// Forward declarations
struct ReproductiveInstruction;
using SpeciesInstructionSet = std::vector<ReproductiveInstruction>;
using GenerationInstructionSets = std::unordered_map<uint32_t, SpeciesInstructionSet>;

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
        std::multimap<FitnessResultType, DynamicGenomeData>& orderedGenomeData,
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const GenerationInstructionSets& instructionSets,
        const DynamicDataUpdateParams& params
    );

    const uint32_t _maxProtectionLimit;           // Max protection counter before elimination
    const uint32_t _maxSpeciesProtectionRating;  // Max species protection rating before elimination
    const double _protectedTierPercentage;       // Percentage of population in protected tier (e.g., 0.3 for 30%)
    const uint32_t _worstSpeciesCount;           // Number of worst species to penalize each generation
    const uint32_t _minActiveSpeciesCount;       // Minimum number of active species before elimination/protection penalties are disabled
};

// Helper function declarations
void updateInstructionSetSizes(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationInstructionSets& instructionSets
);

void updateSpeciesRanking(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, double>& speciesAverageRanks
);

// Main operator function - templated on fitness result type
template<typename FitnessResultType>
void dynamicDataUpdate(
    std::multimap<FitnessResultType, DynamicGenomeData>& orderedGenomeData,  // Best to worst fitness
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationInstructionSets& instructionSets,                       // NEW: Instruction sets from GenerationPlanner
    const DynamicDataUpdateParams& params
) {
    // Phase 1: Instruction Set Counting and Size Updates (NEW)
    updateInstructionSetSizes(speciesData, instructionSets);
    
    // Phase 2: Species Performance Analysis - Single pass with running averages
    std::unordered_map<uint32_t, double> speciesRankSum;
    std::unordered_map<uint32_t, size_t> speciesCount;
    std::unordered_map<uint32_t, double> speciesAverageRanks;
    
    // Track newly discovered species to set their instructionSetsSize after the loop
    std::vector<uint32_t> newlyDiscoveredSpecies;
    
    // Reset species population sizes at start of analysis
    for (auto& [speciesId, data] : speciesData) {
        data.currentPopulationSize = 0;
        speciesRankSum[speciesId] = 0.0;
        speciesCount[speciesId] = 0;
    }
    
    size_t rank = 0;
    for (auto& [fitnessResult, genomeData] : orderedGenomeData) {
        const uint32_t speciesId = genomeData.speciesId;
        
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
            newSpecies.instructionSetsSize = 0;       // Will be set to currentPopulationSize after the loop
            newSpecies.protectionRating = 0;          // Default protection rating
            newSpecies.speciesRank = 0;               // Will be calculated in ranking phase
            newSpecies.isMarkedForElimination = false; // Default elimination status
            speciesData[speciesId] = newSpecies;
            
            // Track this as a newly discovered species
            newlyDiscoveredSpecies.push_back(speciesId);
            // LOG_DEBUG("Created species {}: currentPopulationSize=1, instructionSetsSize=0 (will be set to population size)", speciesId);
        }
        
        ++rank;
    }
    
    // Set instructionSetsSize for newly discovered species (they were copied as-is)
    for (uint32_t speciesId : newlyDiscoveredSpecies) {
        auto& speciesInfo = speciesData[speciesId];
        speciesInfo.instructionSetsSize = speciesInfo.currentPopulationSize;
        // auto logger = LOGGER("population.DynamicDataUpdate");
        // LOG_DEBUG("Set newly discovered species {} instructionSetsSize to {} (equals currentPopulationSize)", 
            // speciesId, speciesInfo.instructionSetsSize);
    }
    
    // Calculate average rank per species
    for (const auto& [speciesId, count] : speciesCount) {
        if (count > 0) {
            speciesAverageRanks[speciesId] = speciesRankSum[speciesId] / count;
        }
    }
    
    // Update species rankings based on performance (NEW)
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
    
    // Phase 3: Individual Genome Protection Updates
    
    // Count active species before applying protection penalties
    uint32_t activeSpeciesCount = 0;
    for (const auto& [speciesId, data] : speciesData) {
        if (!data.isMarkedForElimination) {
            activeSpeciesCount++;
        }
    }
    
    // Only apply protection penalties if we have more than minimum active species
    const bool applyProtectionPenalties = activeSpeciesCount > params._minActiveSpeciesCount;
    
    const size_t totalPopulation = orderedGenomeData.size();
    const size_t protectedTierThreshold = static_cast<size_t>(
        totalPopulation * params._protectedTierPercentage
    );
    
    rank = 0;
    for (auto& [fitnessResult, genomeData] : orderedGenomeData) {
        // Skip protection updates for genomes under repair
        if (genomeData.isUnderRepair) {
            ++rank;
            continue;
        }
        
        // Check if genome is in protected tier (bottom X% by rank)
        const bool isInProtectedTier = rank >= (totalPopulation - protectedTierThreshold);
        
        if (applyProtectionPenalties) {
            if (isInProtectedTier) {
                genomeData.protectionCounter++;
            } else {
                genomeData.protectionCounter = 0;
            }
            
            // Mark for elimination if protection limit exceeded
            if (genomeData.protectionCounter > params._maxProtectionLimit) {
                genomeData.isMarkedForElimination = true;
            }
        }
        // If we're at/below minimum species count, don't increase protection counters or mark genomes for elimination
        
        ++rank;
    }
    
    // Phase 4: Species Protection Rating Updates
    // Create set of worst species IDs for quick lookup
    std::unordered_set<uint32_t> worstSpeciesIds;
    for (size_t i = speciesRankings.size() - worstSpeciesCount; i < speciesRankings.size(); ++i) {
        worstSpeciesIds.insert(speciesRankings[i].first);
    }
    
    // Update protection ratings: penalize worst species, reset others
    // Only apply species protection penalties if we have more than minimum active species
    // if (applyProtectionPenalties) {
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
    
    // Phase 5: State Consistency Validation (Debug assertions)
#ifndef NDEBUG
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        // Protection counter can exceed maxProtectionLimit for custom scaling logic that punishes worst genomes more severely
        // No assertion needed for protectionCounter upper bound
        assert(speciesData.find(genomeData.speciesId) != speciesData.end());
    }
    
    for (const auto& [speciesId, data] : speciesData) {
        // Protection rating can exceed maxSpeciesProtectionRating for adaptive elimination strategies
        assert(data.speciesRank > 0); // Species rank should be assigned
        
        // Population size should match actual count
        size_t actualCount = 0;
        for (const auto& [fr, gd] : orderedGenomeData) {
            if (gd.speciesId == speciesId) actualCount++;
        }
        assert(data.currentPopulationSize == actualCount);
        
        // Instruction set size validation
        auto instructionIt = instructionSets.find(speciesId);
        if (instructionIt != instructionSets.end()) {
            // Species with instruction sets: size should match actual instruction count
            assert(data.instructionSetsSize == instructionIt->second.size());
        } else {
            // Species not in instruction sets: either newly discovered or previously eliminated species being rediscovered
            // instructionSetsSize should be either 0 (no instruction sets yet) or currentPopulationSize (carry-forward behavior)
            assert(data.instructionSetsSize == 0 || data.instructionSetsSize == data.currentPopulationSize);
        }
    }
#endif
}

} // namespace Population