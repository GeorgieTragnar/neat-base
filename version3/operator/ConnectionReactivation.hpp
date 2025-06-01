#pragma once

#include <vector>
#include <memory>

#include "../data/NodeGene.hpp"
#include "../data/ConnectionGene.hpp"
#include "../data/Genome.hpp"

namespace Operator {

class ConnectionReactivationParams;

Genome connectionReactivation(const Genome& genome, const ConnectionReactivationParams& params);

class ConnectionReactivationParams {
public:
    enum class SelectionStrategy {
        RANDOM,          // Random selection from disabled connections
        OLDEST_FIRST,    // Reactivate connection with lowest innovation number
        NEWEST_FIRST     // Reactivate connection with highest innovation number
    };
    
    ConnectionReactivationParams() = delete;
    ConnectionReactivationParams(SelectionStrategy strategy);

protected:
    friend Genome connectionReactivation(const Genome& genome, const ConnectionReactivationParams& params);
    
    const SelectionStrategy _selectionStrategy;
};

}