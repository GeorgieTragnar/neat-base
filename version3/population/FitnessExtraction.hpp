#pragma once

#include "version3/data/PopulationContainer.hpp"
#include <cstdint>
#include <cassert>

namespace Operator {

// Container-based operator to extract fitness values for given genome indices
template<typename FitnessResultType>
std::pair<FitnessResultType, FitnessResultType> fitnessExtraction(
    const PopulationContainer<FitnessResultType>& container,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& generation
);

// Template implementation
template<typename FitnessResultType>
std::pair<FitnessResultType, FitnessResultType> fitnessExtraction(
    const PopulationContainer<FitnessResultType>& container,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& generation
) {
    const auto& fitnessResults = container.getFitnessResults(generation);
    const auto& genomes = container.getGenomes(generation);
    
    // Validate parent indices
    assert(parentAIndex < genomes.size() && "Parent A index out of bounds");
    assert(parentBIndex < genomes.size() && "Parent B index out of bounds");
    
    // Extract fitness values from multimap
    FitnessResultType fitnessA{}, fitnessB{};
    bool foundA = false, foundB = false;
    
    for (const auto& [fitness, globalIndex] : fitnessResults) {
        if (globalIndex == parentAIndex) {
            fitnessA = fitness;
            foundA = true;
        }
        if (globalIndex == parentBIndex) {
            fitnessB = fitness;
            foundB = true;
        }
        if (foundA && foundB) break;
    }
    
    assert(foundA && "Parent A fitness not found in results");
    assert(foundB && "Parent B fitness not found in results");
    
    return std::make_pair(fitnessA, fitnessB);
}

}