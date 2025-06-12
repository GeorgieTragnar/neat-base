#pragma once

#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cassert>

namespace Population {

// Forward declarations
class DynamicDataUpdateParams;

// Dynamic data structures for genome and species tracking
struct DynamicGenomeData {
    uint32_t speciesId;
    uint32_t protectionCounter = 0;
    bool isUnderRepair = false;
    bool isMarkedForElimination = false;
};

struct DynamicSpeciesData {
    uint32_t protectionRating = 0;
    uint32_t currentPopulationSize = 0;
    bool isMarkedForElimination = false;
};

class DynamicDataUpdateParams {
public:
    DynamicDataUpdateParams() = delete;
    DynamicDataUpdateParams(
        uint32_t maxProtectionLimit,
        uint32_t maxSpeciesProtectionRating,
        double protectedTierPercentage,
        uint32_t worstSpeciesCount = 1
    ) : _maxProtectionLimit(maxProtectionLimit),
        _maxSpeciesProtectionRating(maxSpeciesProtectionRating),
        _protectedTierPercentage(protectedTierPercentage),
        _worstSpeciesCount(worstSpeciesCount) {
        
        assert(protectedTierPercentage >= 0.0 && protectedTierPercentage <= 1.0);
        assert(worstSpeciesCount > 0);
    }

protected:
    template<typename FitnessResultType>
    friend void dynamicDataUpdate(
        std::map<FitnessResultType, DynamicGenomeData>& orderedGenomeData,
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const EvolutionLoopFinishingOperatorParams& params
    );

    const uint32_t _maxProtectionLimit;           // Max protection counter before elimination
    const uint32_t _maxSpeciesProtectionRating;  // Max species protection rating before elimination
    const double _protectedTierPercentage;       // Percentage of population in protected tier (e.g., 0.3 for 30%)
    const uint32_t _worstSpeciesCount;           // Number of worst species to penalize each generation
};

// Main operator function - templated on fitness result type
template<typename FitnessResultType>
void dynamicDataUpdate(
    std::map<FitnessResultType, DynamicGenomeData>& orderedGenomeData,  // Best to worst fitness
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EvolutionLoopFinishingOperatorParams& params
) {
    // 1. Species Performance Analysis
    std::unordered_map<uint32_t, std::vector<size_t>> speciesMembers; // speciesId -> ranks
    std::unordered_map<uint32_t, double> speciesAverageRanks;
    
    size_t rank = 0;
    for (auto& [fitnessResult, genomeData] : orderedGenomeData) {
        speciesMembers[genomeData.speciesId].push_back(rank);
        ++rank;
    }
    
    // Calculate average rank per species
    for (const auto& [speciesId, ranks] : speciesMembers) {
        double sum = 0.0;
        for (size_t r : ranks) {
            sum += static_cast<double>(r);
        }
        speciesAverageRanks[speciesId] = sum / ranks.size();
    }
    
    // Rank species from best to worst performing
    std::vector<std::pair<uint32_t, double>> speciesRankings;
    for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
        speciesRankings.emplace_back(speciesId, avgRank);
    }
    std::sort(speciesRankings.begin(), speciesRankings.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // 2. Individual Genome Protection Updates
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
        const bool isInProtectedTier = (totalPopulation - rank) <= protectedTierThreshold;
        
        if (isInProtectedTier) {
            genomeData.protectionCounter++;
        } else {
            genomeData.protectionCounter = 0;
        }
        
        // Mark for elimination if protection limit exceeded
        if (genomeData.protectionCounter > params._maxProtectionLimit) {
            genomeData.isMarkedForElimination = true;
        }
        
        ++rank;
    }
    
    // 3. Species Protection Rating Updates
    const size_t worstSpeciesCount = std::min(
        static_cast<size_t>(params._worstSpeciesCount), 
        speciesRankings.size()
    );
    
    // Create set of worst species IDs for quick lookup
    std::unordered_set<uint32_t> worstSpeciesIds;
    for (size_t i = speciesRankings.size() - worstSpeciesCount; i < speciesRankings.size(); ++i) {
        worstSpeciesIds.insert(speciesRankings[i].first);
    }
    
    // Update protection ratings: penalize worst species, reset others
    for (auto& [speciesId, data] : speciesData) {
        if (worstSpeciesIds.find(speciesId) != worstSpeciesIds.end()) {
            // This species is in worst N - apply penalty
            data.protectionRating++;
            
            // Mark species for elimination if rating threshold exceeded
            if (data.protectionRating > params._maxSpeciesProtectionRating) {
                data.isMarkedForElimination = true;
            }
        } else {
            // This species escaped worst N - reset rating
            data.protectionRating = 0;
        }
    }
    
    // 4. Update Species Population Sizes
    for (auto& [speciesId, data] : speciesData) {
        data.currentPopulationSize = 0;
    }
    
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        if (speciesData.find(genomeData.speciesId) != speciesData.end()) {
            speciesData[genomeData.speciesId].currentPopulationSize++;
        }
    }
    
    // 5. State Consistency Validation (Debug assertions)
    #ifdef DEBUG
    for (const auto& [fitnessResult, genomeData] : orderedGenomeData) {
        assert(genomeData.protectionCounter <= params._maxProtectionLimit + 1); // +1 for the elimination frame
        assert(speciesData.find(genomeData.speciesId) != speciesData.end());
    }
    
    for (const auto& [speciesId, data] : speciesData) {
        assert(data.protectionRating <= params._maxSpeciesProtectionRating + 1); // +1 for the elimination frame
        // Population size should match actual count
        size_t actualCount = 0;
        for (const auto& [fr, gd] : orderedGenomeData) {
            if (gd.speciesId == speciesId) actualCount++;
        }
        assert(data.currentPopulationSize == actualCount);
    }
    #endif
}

}