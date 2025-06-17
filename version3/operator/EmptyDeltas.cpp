#include "EmptyDeltas.hpp"

namespace Operator {

bool emptyDeltas(const Genome& genome) {
    // Access protected delta vectors through friend access
    // Returns true if both delta vectors are empty (clean state)
    // Returns false if either has pending deltas (needs phenotype update)
    return const_cast<Genome&>(genome).get_connectionGeneDeltas().empty() && 
           const_cast<Genome&>(genome).get_nodeGeneDeltas().empty();
}

} // namespace Operator