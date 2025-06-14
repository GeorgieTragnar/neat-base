#pragma once

#include <vector>
#include <cassert>
#include "ReproductiveInstruction.hpp"

namespace Population {

/**
 * @brief Pure data transformation function: resolve relative indices to global indices
 * 
 * Single Responsibility: Transform relative indices to global indices in instruction objects.
 * 
 * Minimal Scope:
 * - Input: Instructions with relative indices + global index lookup table
 * - Operation: Replace relative indices with corresponding global indices  
 * - Output: Same instructions with populated global indices
 * - NO business logic, NO validation, NO context awareness
 * 
 * @param instructions Instructions to modify in-place
 * @param globalIndices Index lookup table (relative index -> global index)
 */
void resolveParentIndices(
    SpeciesInstructionSet& instructions,        // Instructions to modify in-place
    const std::vector<size_t>& globalIndices   // Index lookup table
);

} // namespace Population