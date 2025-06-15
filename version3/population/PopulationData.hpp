#pragma once

#include <cstdint>

class Genome;

namespace Population {

// Dynamic data structures for genome and species tracking
// Shared across population operators for consistent data management

struct DynamicGenomeData {
    uint32_t speciesId = UINT32_MAX;
    uint32_t protectionCounter = 0;
    bool isUnderRepair = false;
    bool isMarkedForElimination = false;
    uint32_t genomeIndex = UINT32_MAX;
};

struct DynamicSpeciesData {
    uint32_t protectionRating = 0;
    uint32_t currentPopulationSize = 0;
    uint32_t instructionSetsSize = 0;  // Number of instruction sets this species will generate (renamed from currentPopulationSize)
    uint32_t speciesRank = 0;          // Ordinal ranking (1st, 2nd, 3rd...) from DynamicDataUpdate - needs implementation
    bool isMarkedForElimination = false;
};

}