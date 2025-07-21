#pragma once

#include "version3/data/PopulationData.hpp"
#include <cstdint>

namespace Operator {

// Creates DynamicGenomeData for crossover offspring
DynamicGenomeData createCrossoverDynamicData(
    const DynamicGenomeData& parentAData,
    const DynamicGenomeData& parentBData,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& speciesId,
    const bool& isUnderRepair
);

// Implementation
inline DynamicGenomeData createCrossoverDynamicData(
    const DynamicGenomeData& parentAData,
    const DynamicGenomeData& parentBData,
    const size_t& parentAIndex,
    const size_t& parentBIndex,
    const uint32_t& speciesId,
    const bool& isUnderRepair
) {
    DynamicGenomeData offspringData;
    offspringData.speciesId = speciesId;  // Use provided species ID (from compatibilityDistance)
    offspringData.pendingEliminationCounter = 0;      // Fresh start
    offspringData.isUnderRepair = isUnderRepair;
    offspringData.isMarkedForElimination = false;
    offspringData.parentAIndex = parentAIndex;
    offspringData.parentBIndex = parentBIndex;
    
    return offspringData;
}

}