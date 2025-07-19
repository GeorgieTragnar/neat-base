#include "DynamicDataUpdate.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

// DynamicDataUpdateParams constructor implementation
DynamicDataUpdateParams::DynamicDataUpdateParams(
    uint32_t maxGenomePendingEliminationLimit,
    uint32_t maxSpeciesPendingEliminationRating,
    double speciesElitePlacementProtectionPercentage,
    double speciesPendingEliminationPercentage,
    double genomesPendingEliminationPercentage,
    uint32_t equilibriumSpeciesCount,
    uint32_t targetPopulationSize
) : _maxGenomePendingEliminationLimit(maxGenomePendingEliminationLimit),
    _maxSpeciesPendingEliminationRating(maxSpeciesPendingEliminationRating),
    _speciesElitePlacementProtectionPercentage(speciesElitePlacementProtectionPercentage),
    _speciesPendingEliminationPercentage(speciesPendingEliminationPercentage),
    _genomesPendingEliminationPercentage(genomesPendingEliminationPercentage),
    _equilibriumSpeciesCount(equilibriumSpeciesCount),
    _targetPopulationSize(targetPopulationSize) {
    
    assert(speciesElitePlacementProtectionPercentage >= 0.0 && speciesElitePlacementProtectionPercentage <= 1.0);
    assert(speciesPendingEliminationPercentage >= 0.0 && speciesPendingEliminationPercentage <= 1.0);
    assert(genomesPendingEliminationPercentage >= 0.0 && genomesPendingEliminationPercentage <= 1.0);
    assert(equilibriumSpeciesCount > 0);
    assert(targetPopulationSize > 0);
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