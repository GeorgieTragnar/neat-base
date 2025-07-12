#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "../data/GlobalIndexRegistry.hpp"
#include "../data/PopulationData.hpp"
#include "../logger/Logger.hpp"

namespace Operator {

class PlotElitesParams;

void plotElites(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotElitesParams& params,
    GlobalIndexRegistry& registry,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
);

class PlotElitesParams {
public:
    PlotElitesParams(
        double elitePercentage,
        size_t minimumElitesPerSpecies,
        size_t maximumElitesPerSpecies
    ) : _elitePercentage(elitePercentage),
        _minimumElitesPerSpecies(minimumElitesPerSpecies),
        _maximumElitesPerSpecies(maximumElitesPerSpecies) {
        
        assert(elitePercentage >= 0.0 && elitePercentage <= 1.0);
        assert(minimumElitesPerSpecies > 0);
        assert(maximumElitesPerSpecies >= minimumElitesPerSpecies);
    }

private:
    friend void plotElites(
        const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
        const PlotElitesParams& params,
        GlobalIndexRegistry& registry,
        const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
    );
    
    const double _elitePercentage;
    const size_t _minimumElitesPerSpecies;
    const size_t _maximumElitesPerSpecies;
};

void plotElites(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotElitesParams& params,
    GlobalIndexRegistry& registry,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData) {
    
    assert(!speciesGroupings.empty());

    auto logger = LOGGER("population.PlotElites");
    
    // Clear all existing elite status
    LOG_DEBUG("ELITE_TRACK: plotElites - Clearing all existing Elite status");
    registry.clearAllEliteStatus();
    
    for (const auto& [speciesId, genomeIndices] : speciesGroupings) {
        assert(!genomeIndices.empty());
        
        // Skip species marked for elimination
        auto speciesDataIt = speciesData.find(speciesId);
        if (speciesDataIt != speciesData.end() && speciesDataIt->second.isMarkedForElimination) {
            continue; // Skip this entire species
        }
        
        // Filter out genomes marked for elimination through global index registry
        std::vector<size_t> activeGenomes;
        for (size_t globalIndex : genomeIndices) {
            if (registry.getState(static_cast<uint32_t>(globalIndex)) == GenomeState::Active) {
                activeGenomes.push_back(globalIndex);
            }
        }
        
        // Note: SpeciesGrouping should not produce duplicates (validated by assertion in SpeciesGrouping.hpp:65)
        // However, we preserve duplicate removal using fitness-order-preserving approach
        std::vector<size_t> uniqueActiveGenomes;
        std::unordered_set<size_t> seenIndices;
        for (size_t globalIndex : activeGenomes) {
            if (seenIndices.find(globalIndex) == seenIndices.end()) {
                uniqueActiveGenomes.push_back(globalIndex);
                seenIndices.insert(globalIndex);
            }
        }
        activeGenomes = std::move(uniqueActiveGenomes);
        
        // Skip species with no active genomes
        if (activeGenomes.empty()) {
            continue;
        }
        
        // Calculate elite count for this species based on active genomes
        size_t speciesSize = activeGenomes.size();
        size_t eliteCount = static_cast<size_t>(std::ceil(speciesSize * params._elitePercentage));
        
        // Apply minimum and maximum constraints
        eliteCount = std::max(eliteCount, params._minimumElitesPerSpecies);
        eliteCount = std::min(eliteCount, params._maximumElitesPerSpecies);
        eliteCount = std::min(eliteCount, speciesSize); // Can't select more elites than active genomes in species
        
        // Mark top performers as elite (active genomes are fitness-ordered from worst to best, so take from end)
        LOG_DEBUG("ELITE_TRACK: plotElites - Species {} selecting {} elites from {} active genomes", 
                 speciesId, eliteCount, speciesSize);
        for (size_t i = 0; i < eliteCount; ++i) {
            uint32_t eliteIndex = static_cast<uint32_t>(activeGenomes[activeGenomes.size() - 1 - i]);
            registry.markAsElite(eliteIndex);
            LOG_DEBUG("ELITE_TRACK: plotElites - Marked genome {} as Elite for species {}", eliteIndex, speciesId);
        }
    }
}

}