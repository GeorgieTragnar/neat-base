#pragma once

#include "../data/Genome.hpp"

namespace Operator {

/**
 * @brief Check if genome has any active (enabled) connections
 * 
 * This operator checks if the genome contains any connection genes
 * that are enabled (enabled = true). This is useful for determining
 * whether node mutation can be applied (requires at least one active connection to split).
 * 
 * @param genome The genome to check
 * @return true if at least one connection is enabled, false if all connections are disabled
 */
bool hasActiveConnections(const Genome& genome);

} // namespace Operator