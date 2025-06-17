#pragma once

#include "../data/Genome.hpp"

namespace Operator {

/**
 * @brief Check if genome has empty deltas (no uncleared mutations)
 * 
 * This operator checks if both connection and node gene deltas are empty,
 * which indicates the genome is in a clean state with no pending mutations
 * that need phenotype updates.
 * 
 * @param genome The genome to check
 * @return true if both connection and node deltas are empty, false otherwise
 */
bool emptyDeltas(const Genome& genome);

} // namespace Operator