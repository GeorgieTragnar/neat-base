#include "ResolveParentIndices.hpp"

namespace Population {

void resolveParentIndices(
    SpeciesInstructionSet& instructions,
    const std::vector<size_t>& globalIndices
) {
    // Simple lookup loop - pure data transformation
    for (auto& instruction : instructions) {
        // Clear any existing global indices
        instruction.globalParentIndices.clear();
        
        // Reserve space for efficiency (avoid reallocations)
        instruction.globalParentIndices.reserve(instruction.relativeParentIndices.size());
        
        // Transform each relative index to global index
        for (uint32_t relativeIndex : instruction.relativeParentIndices) {
            // Debug-only bounds checking (upstream operators responsible for correctness)
            #ifdef DEBUG
            assert(relativeIndex < globalIndices.size() && "Relative index out of bounds");
            #endif
            
            // Direct array access: O(1) per lookup
            instruction.globalParentIndices.push_back(
                static_cast<uint32_t>(globalIndices[relativeIndex])
            );
        }
    }
}

} // namespace Population