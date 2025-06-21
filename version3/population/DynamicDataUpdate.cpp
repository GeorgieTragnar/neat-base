#include "DynamicDataUpdate.hpp"
#include "../logger/Logger.hpp"

namespace Population {

// DynamicDataUpdateParams constructor implementation
DynamicDataUpdateParams::DynamicDataUpdateParams(
    uint32_t maxProtectionLimit,
    uint32_t maxSpeciesProtectionRating,
    double protectedTierPercentage,
    uint32_t worstSpeciesCount,
    uint32_t minActiveSpeciesCount
) : _maxProtectionLimit(maxProtectionLimit),
    _maxSpeciesProtectionRating(maxSpeciesProtectionRating),
    _protectedTierPercentage(protectedTierPercentage),
    _worstSpeciesCount(worstSpeciesCount),
    _minActiveSpeciesCount(minActiveSpeciesCount) {
    
    assert(protectedTierPercentage >= 0.0 && protectedTierPercentage <= 1.0);
    assert(worstSpeciesCount > 0);
    assert(minActiveSpeciesCount > 0);
}


// Internal helper function for species ranking based on average performance
void updateSpeciesRanking(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, double>& speciesAverageRanks
) {
    // Create vector of (speciesId, averageRank) pairs for sorting
    std::vector<std::pair<uint32_t, double>> speciesRankings;
    for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
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
    
    // Handle species not in rankings (no genomes) - assign worst rank
    uint32_t worstRank = static_cast<uint32_t>(speciesRankings.size() + 1);
    for (auto& [speciesId, data] : speciesData) {
        if (data.speciesRank == 0) { // Not yet assigned
            data.speciesRank = worstRank;
        }
    }
}

} // namespace Population