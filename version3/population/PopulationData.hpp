#pragma once

#include <cstdint>

namespace Population {

// Dynamic data structures for genome and species tracking
// Shared across population operators for consistent data management

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

}