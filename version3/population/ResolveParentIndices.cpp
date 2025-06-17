#include "ResolveParentIndices.hpp"
// #include "../logger/Logger.hpp"

namespace Population {

void resolveParentIndices(
    SpeciesInstructionSet& instructions,
    const std::vector<size_t>& globalIndices
) {
    // static auto logger = LOGGER("population.ResolveParentIndices");
    // LOG_DEBUG("resolveParentIndices: Processing {} instructions with {} global indices", 
    //     instructions.size(), globalIndices.size());
    
    // // Log global indices for debugging
    // for (size_t i = 0; i < globalIndices.size(); ++i) {
    //     LOG_DEBUG("  globalIndices[{}] = {}", i, globalIndices[i]);
    // }
    
    // Simple lookup loop - pure data transformation
    for (auto& instruction : instructions) {
        // LOG_DEBUG("Processing instruction with {} relative parent indices", 
        //     instruction.relativeParentIndices.size());
        
        // Clear any existing global indices
        instruction.globalParentIndices.clear();
        
        // Reserve space for efficiency (avoid reallocations)
        instruction.globalParentIndices.reserve(instruction.relativeParentIndices.size());
        
        // Transform each relative index to global index
        for (uint32_t relativeIndex : instruction.relativeParentIndices) {
            // LOG_DEBUG("  Converting relative index {} to global index", relativeIndex);
            
            // Use modulo to wrap around when relative index exceeds actual population size
            size_t actualIndex = relativeIndex % globalIndices.size();
            // LOG_DEBUG("    Relative {} -> Actual {} (modulo {})", relativeIndex, actualIndex, globalIndices.size());
            
            // Direct array access: O(1) per lookup
            uint32_t globalIndex = static_cast<uint32_t>(globalIndices[actualIndex]);
            // LOG_DEBUG("    Final: Relative {} -> Global {}", relativeIndex, globalIndex);
            instruction.globalParentIndices.push_back(globalIndex);
        }
    }
}

} // namespace Population