#pragma once

#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cassert>
#include "../data/PopulationData.hpp"
#include "../data/GlobalIndexRegistry.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

class DynamicDataUpdateParams;

template<typename FitnessResultType>
void dynamicDataUpdate(
    const std::multimap<FitnessResultType, size_t>& fitnessResults,
    std::vector<DynamicGenomeData>& genomeData,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping,
    const DynamicDataUpdateParams& params,
    GlobalIndexRegistry& registry
);

// Helper function declarations
void updateSpeciesRanking(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, double>& speciesAverageRanks
);

class DynamicDataUpdateParams {
public:
    DynamicDataUpdateParams() = delete;
    DynamicDataUpdateParams(
        uint32_t maxGenomePendingEliminationLimit,
        uint32_t maxSpeciesPendingEliminationRating,
        double speciesElitePlacementProtectionPercentage,
        double speciesPendingEliminationPercentage,
        double genomesPendingEliminationPercentage,
        uint32_t equilibriumSpeciesCount,
        uint32_t targetPopulationSize
    );

protected:
    template<typename FitnessResultType>
    friend void dynamicDataUpdate(const std::multimap<FitnessResultType, size_t>& fitnessResults, std::vector<DynamicGenomeData>& genomeData, std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData, const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping, const DynamicDataUpdateParams& params, GlobalIndexRegistry& registry);

    const uint32_t _maxGenomePendingEliminationLimit;           // Max pending elimination counter before actual elimination
    const uint32_t _maxSpeciesPendingEliminationRating;  // Max species pending elimination rating before elimination
    const double _speciesElitePlacementProtectionPercentage; // special cutoff protection for species if their champion is in a certain percentage of the population
    const double _speciesPendingEliminationPercentage; // what percentage of excess species above equilibrium species count is susceptible to pending elimination
    const double _genomesPendingEliminationPercentage; // what percentage of excess genomes inside single species is above equilibrium value (targetPopulation / activespeciescount) is susceptible to pending elimination
    const uint32_t _equilibriumSpeciesCount;               // Minimum number of active species before elimination penalties are disabled
    const uint32_t _targetPopulationSize;                // Target population size for equilibrium calculation
};

// Template function implementation
template<typename FitnessResultType>
void dynamicDataUpdate(
    const std::multimap<FitnessResultType, size_t>& fitnessResults,  // Best to worst fitness with global indices
    std::vector<DynamicGenomeData>& genomeData,  // Genome data vector for direct access
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping,  // Species grouping data
    const DynamicDataUpdateParams& params,
    GlobalIndexRegistry& registry
) {
    // Phase 1: Species Performance Analysis - Single pass with running averages
    auto logger = LOGGER("population.DynamicDataUpdate");
    // LOG_DEBUG("ELIMINATION PHASE 1: Starting species performance analysis for {} genomes", fitnessResults.size());
    
    std::unordered_map<uint32_t, size_t> speciesRankSum;
    std::unordered_map<uint32_t, size_t> speciesCount;
    std::unordered_map<uint32_t, double> speciesAverageRanks;
    std::unordered_map<uint32_t, size_t> speciesBestGenomeRank;
    
    // Reset species population sizes at start of analysis
    // LOG_DEBUG("ELIMINATION PHASE 1: Resetting {} species population counters", speciesData.size());
    size_t speciesWithEliminated = 0;
    for (auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination) speciesWithEliminated++;
        data.currentPopulationSize = 0;
    }
    // LOG_DEBUG("ELIMINATION PHASE 1: Found {} species currently marked for elimination", speciesWithEliminated);
    
    size_t rank = 0;
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        
        // Note: Species discovery happens below - genome can reference previously unknown species
        // This is valid when new species emerge during evolution
        
        // Update species rank sum and count
        speciesRankSum[speciesId] += rank;
        speciesCount[speciesId]++;
        
        // Track best genome rank for each species (fitness results are ordered worst to best)
        // Since we iterate worst to best, the last rank seen for each species is their best performer
        speciesBestGenomeRank[speciesId] = rank;
        
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
            newSpecies.pendingEliminationRating = 0;          // Default protection rating
            newSpecies.speciesRank = 0;               // TODO: Review new species rank assignment timing
            newSpecies.isMarkedForElimination = false; // Default elimination status
            speciesData[speciesId] = newSpecies;
        }
        
        ++rank;
    }
    
    // Calculate average rank per species
    // LOG_DEBUG("ELIMINATION PHASE 1: Calculating average ranks for {} species", speciesCount.size());
    for (const auto& [speciesId, count] : speciesCount) {
        if (count > 0) {
            speciesAverageRanks[speciesId] = static_cast<double>(speciesRankSum[speciesId]) / count;
            // LOG_TRACE("Species {}: {} genomes, avg rank {:.2f}", speciesId, count, speciesAverageRanks[speciesId]);
        }
    }

    uint32_t activeSpeciesCount = 0;
    for (const auto& [speciesId, data] : speciesData) {
        if (data.currentPopulationSize > 0)
            ++activeSpeciesCount;
    }
    
#ifndef NDEBUG
    // Log species performance summary
    std::vector<std::pair<uint32_t, double>> speciesPerf;
    for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
        speciesPerf.emplace_back(speciesId, avgRank);
    }
    std::sort(speciesPerf.begin(), speciesPerf.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    
    std::string speciesPerfStr = "[";
    for (size_t i = 0; i < speciesPerf.size(); ++i) {
        if (i > 0) speciesPerfStr += ", ";
        speciesPerfStr += fmt::format("{}({:.1f})", speciesPerf[i].first, speciesPerf[i].second);
    }
    speciesPerfStr += "]";
    LOG_DEBUG("ELIMINATION PHASE 1: Species performance (best to worst): {}", speciesPerfStr);
#endif
    
    // Update species rankings based on performance
    updateSpeciesRanking(speciesData, speciesAverageRanks);

#ifndef NDEBUG
    for (auto& [speciesId, data] : speciesData) {
        assert(data.speciesRank > 0 && "species rank must be positive and non zero");
    }
#endif
    
// Phase 2: Individual Genome Pending Elimination Updates
    // LOG_DEBUG("ELIMINATION PHASE 2: Starting individual genome pending elimination updates");
    
    
    const uint32_t excessSpeciesCount = (activeSpeciesCount > params._equilibriumSpeciesCount)
        ? activeSpeciesCount - params._equilibriumSpeciesCount : 0;
    const uint32_t equilibrium = (activeSpeciesCount == 0) ? 0 :
        (activeSpeciesCount > params._targetPopulationSize) ? 1 : params._targetPopulationSize / activeSpeciesCount;

    assert(excessSpeciesCount < activeSpeciesCount && "there cannot be more excess than active species");

    const uint32_t speciesPendingElimination = excessSpeciesCount > 0
        ? excessSpeciesCount * params._speciesPendingEliminationPercentage : 0;

    assert(speciesPendingElimination < activeSpeciesCount && "there cannot be more species pending elimination that active species");

    // Species elimination logic: Check if worst performing species need pending elimination
    std::unordered_set<uint32_t> speciesWithIncreasedRating;
    std::unordered_set<uint32_t> eliminatedSpeciesIds;
    
    if (speciesPendingElimination > 0) {
        // Create vector of species sorted by average rank (worst to best performance)
        std::vector<std::pair<uint32_t, double>> speciesRanking;
        for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
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
        const size_t totalPopulation = fitnessResults.size();
        const size_t speciesToCheck = std::min(static_cast<size_t>(speciesPendingElimination), 
                                               speciesRanking.size());
        
        for (size_t i = 0; i < speciesToCheck; ++i) {
            const uint32_t speciesId = speciesRanking[i].first;
            
            // Get the best genome rank for this species (their champion)
            auto championRankIt = speciesBestGenomeRank.find(speciesId);
            if (championRankIt == speciesBestGenomeRank.end()) {
                // This shouldn't happen, but handle gracefully
                LOG_WARN("Species {} has no champion rank recorded, skipping elimination check", speciesId);
                continue;
            }
            
            const size_t championRank = championRankIt->second;
            
            // Calculate champion's percentile in population (0.0 = worst, 1.0 = best)
            const double championPercentile = (totalPopulation > 1) ? 
                static_cast<double>(championRank) / (totalPopulation - 1) : 1.0;
            
            // Check if champion is protected by elite placement percentage
            if (championPercentile >= params._speciesElitePlacementProtectionPercentage) {
                // Champion is in top X% of population - species is protected from elimination
                LOG_DEBUG("Species {} protected: champion at rank {} (percentile {:.1f}%) >= protection threshold {:.1f}%", 
                         speciesId, championRank, championPercentile * 100.0, 
                         params._speciesElitePlacementProtectionPercentage * 100.0);
                continue;
            }
            
            // Species is not protected - increase pending elimination rating
            auto& speciesInfo = speciesData[speciesId];
            speciesInfo.pendingEliminationRating++;
            speciesWithIncreasedRating.insert(speciesId);
            
            LOG_DEBUG("Species {} pending elimination increased to {}: champion at rank {} (percentile {:.1f}%) < protection threshold {:.1f}%", 
                     speciesId, speciesInfo.pendingEliminationRating, championRank, 
                     championPercentile * 100.0, params._speciesElitePlacementProtectionPercentage * 100.0);
            
            // Mark species for elimination if rating exceeds limit
            if (speciesInfo.pendingEliminationRating >= params._maxSpeciesPendingEliminationRating) {
                // ELITE_TRACK: Check if species being eliminated contains Elite genomes
                bool hasEliteGenomes = false;
                for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
                    if (genomeData[globalIndex].speciesId == speciesId && 
                        registry.getState(globalIndex) == GenomeState::Elite) {
                        hasEliteGenomes = true;
                        LOG_DEBUG("ELITE_TRACK: DynamicDataUpdate - Species {} has Elite genome {} but being eliminated!", 
                                 speciesId, globalIndex);
                    }
                }
                
                speciesInfo.isMarkedForElimination = true;
                eliminatedSpeciesIds.insert(speciesId);
                LOG_INFO("Species {} marked for elimination: pending rating {} >= limit {} (hasElites: {})", 
                        speciesId, speciesInfo.pendingEliminationRating, params._maxSpeciesPendingEliminationRating, hasEliteGenomes);
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
    
    // Mark all genomes belonging to eliminated species for elimination
    if (!eliminatedSpeciesIds.empty()) {
        size_t genomesMarkedForElimination = 0;
        for (size_t globalIndex = 0; globalIndex < genomeData.size(); ++globalIndex) {
            auto& genome = genomeData[globalIndex];
            if (eliminatedSpeciesIds.find(genome.speciesId) != eliminatedSpeciesIds.end()) {
                if (registry.getState(globalIndex) == GenomeState::Active) {
                    registry.markForElimination(globalIndex);
                    genomesMarkedForElimination++;
                }
            }
        }
        LOG_INFO("Marked {} genomes for elimination across {} eliminated species", 
                genomesMarkedForElimination, eliminatedSpeciesIds.size());
    }
    
    // for each active species identify genomes pending elimination
    for (auto& [speciesId, data] : speciesData) {
        // extinct, marked for elimination or without valid genomes species are skipped
        if (data.currentPopulationSize == 0 || data.isMarkedForElimination 
            || speciesGrouping.find(speciesId) == speciesGrouping.end())
            continue;

        const std::vector<size_t>& validGenomes = speciesGrouping.at(speciesId);
        const uint32_t activeGenomeCount = static_cast<uint32_t>(validGenomes.size());
        const uint32_t excessGenomeCount = (activeGenomeCount > equilibrium)
            ? activeGenomeCount - equilibrium : 0;

        bool skip = excessGenomeCount == 0;

        assert(excessGenomeCount < activeGenomeCount && "we cannot have more excess genomes that the active species population size");

        const uint32_t genomesPendingElimination = excessGenomeCount > 0
            ? excessGenomeCount * params._genomesPendingEliminationPercentage : 0;

        skip = skip || genomesPendingElimination == 0;

        assert(genomesPendingElimination < activeGenomeCount && "we cannot mark more genomes for pending elimination than the active species population size");

        assert(!validGenomes.empty() && "at this point species without valid genomes should have been skipped already");

        // increment genomes that are pending elimination
        for (size_t i = validGenomes.size() - genomesPendingElimination;i < validGenomes.size(); ++i) {
            // Skip genomes that may have been marked for elimination during this update cycle
            if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                continue;
            }
            
            auto& data = genomeData[validGenomes[i]];

            assert(!data.isUnderRepair && "invalid access to genome data under repair");

            if (data.pendingEliminationCounter++ >= params._maxGenomePendingEliminationLimit) {
                registry.markForElimination(validGenomes[i]);
            }
        }

        // decrement genomes that are not pending elimination
        for (size_t i = 0; i < (validGenomes.size() - genomesPendingElimination); ++i) {
            // Skip genomes that may have been marked for elimination during this update cycle
            if (registry.getState(validGenomes[i]) != GenomeState::Active) {
                continue;
            }
            
            auto& data = genomeData[validGenomes[i]];

            assert(!data.isUnderRepair && "invalid access to genome data under repair");

            if (data.pendingEliminationCounter != 0)
                data.pendingEliminationCounter--;
        }
    }

#ifndef NDEBUG
    // State Consistency Validation and Logging
    LOG_DEBUG("Starting state consistency validation");
    
    // Enhanced validation with detailed logging
    size_t validationErrors = 0;
    
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        // Validate species existence
        if (speciesData.find(genomeData[globalIndex].speciesId) == speciesData.end()) {
            LOG_ERROR("VALIDATION ERROR: Genome {} references non-existent species {}", 
                     globalIndex, genomeData[globalIndex].speciesId);
            validationErrors++;
        }
        
        // Validate elimination logic consistency - only check genomes newly marked in this cycle
        const auto& genome = genomeData[globalIndex];
        if (registry.getState(globalIndex) == GenomeState::HotElimination && !genome.isUnderRepair) {
            // Check if genome was eliminated via species elimination or individual elimination
            bool eliminatedViaSpecies = eliminatedSpeciesIds.find(genome.speciesId) != eliminatedSpeciesIds.end();
            bool eliminatedViaCounter = genome.pendingEliminationCounter >= params._maxGenomePendingEliminationLimit;
            bool eliminatedViaRepair = genome.repairAttempts > 0; // Could be eliminated by RepairOperator
            
            if (!eliminatedViaSpecies && !eliminatedViaCounter && !eliminatedViaRepair) {
                LOG_ERROR("VALIDATION ERROR: Genome {} marked for elimination but neither species eliminated nor counter {} >= limit {} nor repair attempts {}", 
                         globalIndex, genome.pendingEliminationCounter, params._maxGenomePendingEliminationLimit, genome.repairAttempts);
                validationErrors++;
            }
        }
    }
    
    // Validate species population counts
    for (const auto& [speciesId, data] : speciesData) {
        // Pending elimination rating can exceed maxSpeciesPendingEliminationRating for adaptive elimination strategies
        // Note: speciesRank can be 0 if a species has a single genome that is the best performer (rank 0)
        // Species rank is calculated from average genome ranks, so single-genome species inherit their genome's rank
        if (data.speciesRank < 0) {
            LOG_ERROR("VALIDATION ERROR: Species {} has invalid rank {}", speciesId, data.speciesRank);
            validationErrors++;
        }
        
        // Population size should match actual count
        size_t actualCount = 0;
        for (const auto& [fr, globalIndex] : fitnessResults) {
            if (genomeData[globalIndex].speciesId == speciesId) actualCount++;
        }
        if (data.currentPopulationSize != actualCount) {
            LOG_ERROR("VALIDATION ERROR: Species {} recorded population {} != actual count {}", 
                     speciesId, data.currentPopulationSize, actualCount);
            validationErrors++;
        }
    }
    
    LOG_DEBUG("ELIMINATION PHASE 4 COMPLETE: Validation errors: {}", validationErrors);
    
    assert(validationErrors == 0 && "State validation failed - check debug logs for details");
#endif
}

} // namespace Operator