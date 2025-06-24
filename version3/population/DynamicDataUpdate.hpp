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
    friend void dynamicDataUpdate(
        const std::multimap<FitnessResultType, size_t>& fitnessResults,
        std::vector<DynamicGenomeData>& genomeData,
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping,
        const DynamicDataUpdateParams& params
    );

    const uint32_t _maxGenomePendingEliminationLimit;           // Max pending elimination counter before actual elimination
    const uint32_t _maxSpeciesPendingEliminationRating;  // Max species pending elimination rating before elimination
    const double _speciesElitePlacementProtectionPercentage; // special cutoff protection for species if their champion is in a certain percentage of the population
    const double _speciesPendingEliminationPercentage; // what percentage of excess species above equilibrium species count is susceptible to pending elimination
    const double _genomesPendingEliminationPercentage; // what percentage of excess genomes inside single species is above equilibrium value (targetPopulation / activespeciescount) is susceptible to pending elimination
    const uint32_t _equilibriumSpeciesCount;               // Minimum number of active species before elimination penalties are disabled
    const uint32_t _targetPopulationSize;                // Target population size for equilibrium calculation
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
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping,  // Species grouping data
    const DynamicDataUpdateParams& params
) {
    // Phase 1: Species Performance Analysis - Single pass with running averages
    auto logger = LOGGER("population.DynamicDataUpdate");
    LOG_DEBUG("ELIMINATION PHASE 1: Starting species performance analysis for {} genomes", fitnessResults.size());
    
    std::unordered_map<uint32_t, size_t> speciesRankSum;
    std::unordered_map<uint32_t, size_t> speciesCount;
    std::unordered_map<uint32_t, double> speciesAverageRanks;
    
    // Reset species population sizes at start of analysis
    LOG_DEBUG("ELIMINATION PHASE 1: Resetting {} species population counters", speciesData.size());
    size_t speciesWithEliminated = 0;
    for (auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination) speciesWithEliminated++;
        data.currentPopulationSize = 0;
    }
    LOG_DEBUG("ELIMINATION PHASE 1: Found {} species currently marked for elimination", speciesWithEliminated);
    
    size_t rank = 0;
    for (const auto& [fitnessResult, globalIndex] : fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        
        // Note: Species discovery happens below - genome can reference previously unknown species
        // This is valid when new species emerge during evolution
        
        // Update species rank sum and count
        speciesRankSum[speciesId] += rank;
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
            newSpecies.pendingEliminationRating = 0;          // Default protection rating
            newSpecies.speciesRank = 0;               // TODO: Review new species rank assignment timing
            newSpecies.isMarkedForElimination = false; // Default elimination status
            speciesData[speciesId] = newSpecies;
        }
        
        ++rank;
    }
    
    // Calculate average rank per species
    LOG_DEBUG("ELIMINATION PHASE 1: Calculating average ranks for {} species", speciesCount.size());
    for (const auto& [speciesId, count] : speciesCount) {
        if (count > 0) {
            speciesAverageRanks[speciesId] = speciesRankSum[speciesId] / count;
            LOG_TRACE("Species {}: {} genomes, avg rank {:.2f}", speciesId, count, speciesAverageRanks[speciesId]);
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
    LOG_DEBUG("ELIMINATION PHASE 2: Starting individual genome pending elimination updates");
    
    
    const uint32_t excessSpeciesCount = (activeSpeciesCount > params._equilibriumSpeciesCount)
        ? activeSpeciesCount - params._equilibriumSpeciesCount : 0;
    const uint32_t equilibrium = (activeSpeciesCount > params._targetPopulationSize)
        ? 1 : params._targetPopulationSize / activeSpeciesCount;

    assert(excessSpeciesCount < activeSpeciesCount && "there cannot be more excess than active species");

    const uint32_t speciesPendingElimination = excessSpeciesCount > 0
        ? excessSpeciesCount * params._speciesPendingEliminationPercentage : 0;

    assert(speciesPendingElimination < activeSpeciesCount && "there cannot be more species pending elimination that active species");

    // TODO: if speciesPendingElimination is higher than 0 then find this amount
    // of species inside speciesAverageRanks where you first order 
    // speciesPendingElimination amount worst species and check if their worst performing elite
    // is inside top speciesElitePlacementProtectionPercentage percentage of the whole population
    // if not then increase the pending elimination ratings
    // this will require careful design alteration and assesment before implementation
    // as we cannot just use the param value from PlotELites operator
    // but instead we need a way to mark genomes as elites inside genome dynamic data through the
    // elite operator - before that is working we cant safely mark species for elimination
    
    // for each active species identify genomes pending elimination
    for (auto& [speciesId, data] : speciesData) {
        // extinct, marked for elimination or without valid genomes species are skipped
        if (data.currentPopulationSize == 0 || data.isMarkedForElimination 
            || speciesGrouping.find(speciesId) == speciesGrouping.end())
            continue;

        const uint32_t excessGenomeCount = (data.currentPopulationSize > equilibrium)
            ? data.currentPopulationSize - equilibrium : 0;

        bool skip = excessGenomeCount == 0;

        assert(excessGenomeCount < data.currentPopulationSize && "we cannot have more excess genomes that the species population size");

        const uint32_t genomesPendingElimination = excessGenomeCount > 0
            ? excessGenomeCount * params._genomesPendingEliminationPercentage : 0;

        skip = skip || genomesPendingElimination == 0;

        assert(genomesPendingElimination < data.currentPopulationSize && "we cannot mark more genomes for pending elimination than the species population size");

        const std::vector<size_t>& validGenomes = speciesGrouping.at(speciesId);

        assert(!validGenomes.empty() && "at this point species without valid genomes should have been skipped already");

        // increment genomes that are pending elimination
        for (size_t i = validGenomes.size() - genomesPendingElimination;i < validGenomes.size(); ++i) {
            auto& data = genomeData[validGenomes[i]];

            assert(!data.isUnderRepair && "invalid access to genome data under repair");
            assert(!data.isMarkedForElimination && "invalid access to genome data marked for elimination");

            if (data.pendingEliminationCounter++ >= params._maxGenomePendingEliminationLimit) {
                data.isMarkedForElimination = true;
            }
        }

        // decrement genomes that are not pending elimination
        for (size_t i = 0; i < (validGenomes.size() - genomesPendingElimination); ++i) {
            auto& data = genomeData[validGenomes[i]];

            assert(!data.isUnderRepair && "invalid access to genome data under repair");
            assert(!data.isMarkedForElimination && "invalid access to genome data marked for elimination");

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
        
        // Validate elimination logic consistency
        const auto& genome = genomeData[globalIndex];
        if (genome.isMarkedForElimination && genome.pendingEliminationCounter <= params._maxGenomePendingEliminationLimit && !genome.isUnderRepair) {
            LOG_ERROR("VALIDATION ERROR: Genome {} marked for elimination but pending counter {} <= limit {}", 
                     globalIndex, genome.pendingEliminationCounter, params._maxGenomePendingEliminationLimit);
            validationErrors++;
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

} // namespace Population