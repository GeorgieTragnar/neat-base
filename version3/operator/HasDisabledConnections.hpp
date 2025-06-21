#pragma once

#include "../data/Genome.hpp"

namespace Operator {

/**
 * @brief Check if genome has any disabled connections
 * 
 * This operator checks if the genome contains any connection genes
 * that are disabled (enabled = false). This is useful for determining
 * whether connection reactivation mutations can be applied.
 * 
 * @param genome The genome to check
 * @return true if at least one connection is disabled, false if all connections are enabled
 */
bool hasDisabledConnections(const Genome& genome);

} // namespace Operator