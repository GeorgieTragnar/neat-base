#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "GlobalIndexRegistry.hpp"
#include "PopulationData.hpp"

namespace Population {

class PlotElitesParams;

std::vector<size_t> plotElites(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotElitesParams& params,
    const GlobalIndexRegistry& registry,
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
    friend std::vector<size_t> plotElites(
        const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
        const PlotElitesParams& params,
        const GlobalIndexRegistry& registry,
        const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData
    );
    
    const double _elitePercentage;
    const size_t _minimumElitesPerSpecies;
    const size_t _maximumElitesPerSpecies;
};

inline std::vector<size_t> plotElites(
    const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGroupings,
    const PlotElitesParams& params,
    const GlobalIndexRegistry& registry,
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData) {
    
    assert(!speciesGroupings.empty());
    
    std::vector<size_t> eliteIndices;
    
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
        
        // Remove duplicates from activeGenomes (in case speciesGrouping contains duplicates)
        std::sort(activeGenomes.begin(), activeGenomes.end());
        activeGenomes.erase(std::unique(activeGenomes.begin(), activeGenomes.end()), activeGenomes.end());
        
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
        
        // Select top performers (active genomes are fitness-ordered from worst to best, so take from end)
        for (size_t i = 0; i < eliteCount; ++i) {
            eliteIndices.push_back(activeGenomes[activeGenomes.size() - 1 - i]);
        }
    }
    
#ifndef NDEBUG
    // Validate no duplicate indices
    std::vector<size_t> sortedElites = eliteIndices;
    std::sort(sortedElites.begin(), sortedElites.end());
    auto uniqueEnd = std::unique(sortedElites.begin(), sortedElites.end());
    assert(uniqueEnd == sortedElites.end() && "Duplicate elite indices detected");
#endif
    
    return eliteIndices;
}

}