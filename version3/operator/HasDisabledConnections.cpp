#include "HasDisabledConnections.hpp"

namespace Operator {

bool hasDisabledConnections(const Genome& genome) {
    // Access protected connection genes through friend access
    // Check if any connection gene has enabled = false
    const auto& connections = const_cast<Genome&>(genome).get_connectionGenes();
    
    for (const auto& conn : connections) {
        if (!conn.get_attributes().enabled) {
            return true;
        }
    }
    
    return false;
}

} // namespace Operator