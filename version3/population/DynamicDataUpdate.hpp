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
        const DynamicDataUpdateParams& params
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
    const DynamicDataUpdateParams& params
) {
    // 1. Species Performance Analysis - Single pass with running averages
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
    for (auto& [fitnessResult, genomeData] : orderedGenomeData) {
        const uint32_t speciesId = genomeData.speciesId;
        
        // Update species rank sum and count
        speciesRankSum[speciesId] += static_cast<double>(rank);
        speciesCount[speciesId]++;
        
        // Update species population size in single pass
        if (speciesData.find(speciesId) != speciesData.end()) {
            speciesData[speciesId].currentPopulationSize++;
        }
        
        ++rank;
    }
    
    // Calculate average rank per species
    for (const auto& [speciesId, count] : speciesCount) {
        if (count > 0) {
            speciesAverageRanks[speciesId] = speciesRankSum[speciesId] / count;
        }
    }
    
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
        const bool isInProtectedTier = rank >= (totalPopulation - protectedTierThreshold);
        
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
    
    // 4. Species population sizes already updated in step 1
    
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