#pragma once

#include <unordered_map>
#include <vector>
#include <map>
#include <cstdint>

template<typename FitnessResultType>
class GenerationPopulationData {
public:
    // Core grouping from SpeciesGrouping operator
    // Maps species ID to vector of global indices (fitness-ordered, worst to best)
    std::unordered_map<uint32_t, std::vector<uint32_t>> speciesGrouping;
    
    // Species performance metrics (filled by SpeciesRanking)
    std::unordered_map<uint32_t, size_t> speciesRankSum;      // Sum of ranks for each species
    std::unordered_map<uint32_t, size_t> speciesCount;        // Number of genomes per species
    std::unordered_map<uint32_t, double> speciesAverageRanks; // Average rank per species
    std::unordered_map<uint32_t, size_t> speciesBestGenomeRank; // Best (highest) rank per species
    
    // Active species count (calculated by SpeciesRanking)  
    uint32_t activeSpeciesCount;
    
    // Reference to fitness results (passed through from PopulationContainer)
    // Multimap ordered worst to best fitness
    const std::multimap<FitnessResultType, uint32_t>& fitnessResults;
    
    // Constructor - takes fitness results reference
    explicit GenerationPopulationData(const std::multimap<FitnessResultType, uint32_t>& results) 
        : fitnessResults(results), activeSpeciesCount(0) {}
    
    // Disable copy constructor and assignment to prevent accidental copies
    GenerationPopulationData(const GenerationPopulationData&) = delete;
    GenerationPopulationData& operator=(const GenerationPopulationData&) = delete;
    
    // Enable move constructor and assignment
    GenerationPopulationData(GenerationPopulationData&&) = default;
    GenerationPopulationData& operator=(GenerationPopulationData&&) = default;
};