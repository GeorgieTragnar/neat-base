#pragma once

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>
#include "../data/PopulationData.hpp"
#include "../data/PopulationContainer.hpp"
#include "GenerationPopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

// SpeciesRanking operator: Calculates species performance metrics and rankings
template<typename FitnessResultType>
void speciesRanking(
    GenerationPopulationData<FitnessResultType>& populationData,
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
);

// Template implementation
template<typename FitnessResultType>
void speciesRanking(
    GenerationPopulationData<FitnessResultType>& populationData,
    const PopulationContainer<FitnessResultType>& container,
    uint32_t generation,
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
) {
    auto logger = LOGGER("population.SpeciesRanking");
    
    // Get genome data from container
    const auto& genomeData = container.getGenomeData(generation);
    
    // Phase 1: Calculate species performance metrics
    size_t rank = 0;
    for (const auto& [fitnessResult, globalIndex] : populationData.fitnessResults) {
        const uint32_t speciesId = genomeData[globalIndex].speciesId;
        
        // Handle species discovery - create species data if it doesn't exist
        if (speciesData.find(speciesId) == speciesData.end()) {
            LOG_DEBUG("Discovered new species {} in genome data, creating species entry", speciesId);
            DynamicSpeciesData newSpecies;
            newSpecies.currentPopulationSize = 0;          // Will be updated by SpeciesGrouping
            newSpecies.pendingEliminationRating = 0;       // Default protection rating
            newSpecies.speciesRank = 0;                    // Will be calculated below
            newSpecies.isMarkedForElimination = false;     // Default elimination status
            speciesData[speciesId] = newSpecies;
        }
        
        // Update species rank sum and count for active genomes only
        if (populationData.speciesGrouping.find(speciesId) != populationData.speciesGrouping.end()) {
            populationData.speciesRankSum[speciesId] += rank;
            populationData.speciesCount[speciesId]++;
            
            // Track best genome rank for each species (fitness results are ordered worst to best)
            // Since we iterate worst to best, the last rank seen for each species is their best performer
            populationData.speciesBestGenomeRank[speciesId] = rank;
        }
        
        ++rank;
    }
    
    // Phase 2: Calculate average rank per species
    for (const auto& [speciesId, count] : populationData.speciesCount) {
        if (count > 0) {
            populationData.speciesAverageRanks[speciesId] = 
                static_cast<double>(populationData.speciesRankSum[speciesId]) / count;
        }
    }
    
    // Phase 3: Count active species and assign ordinal rankings
    populationData.activeSpeciesCount = 0;
    for (const auto& [speciesId, data] : speciesData) {
        if (data.currentPopulationSize > 0) {
            ++populationData.activeSpeciesCount;
        }
    }
    
    // Phase 4: Sort species by average rank and assign ordinal rankings
    // Create vector of (speciesId, averageRank) pairs for sorting
    std::vector<std::pair<uint32_t, double>> speciesRankings;
    for (const auto& [speciesId, avgRank] : populationData.speciesAverageRanks) {
        speciesRankings.emplace_back(speciesId, avgRank);
    }
    
    // Sort by average rank (lower rank = better performance)
    std::sort(speciesRankings.begin(), speciesRankings.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Assign ordinal rankings (1-based)
    for (size_t i = 0; i < speciesRankings.size(); ++i) {
        uint32_t speciesId = speciesRankings[i].first;
        auto speciesIt = speciesData.find(speciesId);
        if (speciesIt != speciesData.end()) {
            speciesIt->second.speciesRank = static_cast<uint32_t>(i + 1); // 1-based ranking
        }
    }
    
    // Handle species not in rankings (no active genomes) - assign worst rank
    uint32_t worstRank = static_cast<uint32_t>(speciesRankings.size() + 1);
    for (auto& [speciesId, data] : speciesData) {
        if (data.speciesRank == 0) { // Not yet assigned
            data.speciesRank = worstRank;
        }
    }

#ifndef NDEBUG
    // Debug logging: species performance summary
    std::vector<std::pair<uint32_t, double>> speciesPerf;
    for (const auto& [speciesId, avgRank] : populationData.speciesAverageRanks) {
        speciesPerf.emplace_back(speciesId, avgRank);
    }
    std::sort(speciesPerf.begin(), speciesPerf.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    std::string speciesPerfStr = "[";
    for (size_t i = 0; i < speciesPerf.size(); ++i) {
        if (i > 0) speciesPerfStr += ", ";
        speciesPerfStr += std::to_string(speciesPerf[i].first) + "(" + 
                         std::to_string(speciesPerf[i].second) + ")";
    }
    speciesPerfStr += "]";
    LOG_DEBUG("Species performance (best to worst): {}", speciesPerfStr);
    
    // Validation: all species should have positive ranks
    for (const auto& [speciesId, data] : speciesData) {
        assert(data.speciesRank > 0 && "species rank must be positive and non zero");
    }
#endif
}

} // namespace Operator