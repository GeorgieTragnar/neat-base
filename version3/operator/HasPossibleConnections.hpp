#pragma once

#include "../data/Genome.hpp"

namespace Operator {

/**
 * @brief Check if genome has any possible new connections that can be added
 * 
 * This operator checks if the genome can accommodate additional connection genes
 * without violating topological constraints or creating duplicate connections.
 * This is useful for determining whether connection mutation can be applied.
 * 
 * The check considers:
 * - Existing connections (no duplicates allowed)
 * - Node type constraints (OUTPUT nodes cannot be sources, INPUT nodes cannot be targets)
 * - Self-connection restrictions (no node can connect to itself)
 * 
 * @param genome The genome to check
 * @return true if at least one new connection is possible, false if genome is fully connected
 */
bool hasPossibleConnections(const Genome& genome);

} // namespace Operator