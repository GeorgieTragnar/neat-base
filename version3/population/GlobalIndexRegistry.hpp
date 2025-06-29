#pragma once

#include <vector>
#include <cstdint>
#include <limits>

namespace Population {

constexpr uint32_t INVALID_INDEX = std::numeric_limits<uint32_t>::max();

enum class GenomeState {
    Active,               // Normal genome that can be used as parent
    HotElimination,       // Marked for elimination, offspring may still exist
    ColdElimination,      // All offspring processed, safe to clean
    ReadyForReplacement   // Cleaned, index can be reused
};

class GlobalIndexRegistry {
public:
    explicit GlobalIndexRegistry(uint32_t maxIndices);

    // State queries
    GenomeState getState(uint32_t globalIndex) const;
    
    // State transitions
    void markForElimination(uint32_t globalIndex);
    void transitionToCold(uint32_t globalIndex);
    void markReadyForReplacement(uint32_t globalIndex);
    void resetToActive(uint32_t globalIndex);
    
    // Index management
    uint32_t getFreeIndex();  // Returns INVALID_INDEX if none available
    void incrementMaxIndex();  // Extend the indexing range
    
    uint32_t getMaxIndex() const { return static_cast<uint32_t>(_states.size()); }

private:
    std::vector<GenomeState> _states;
};

} // namespace Population