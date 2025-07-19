#pragma once

#include <vector>
#include <memory>

#include "version3/data/NodeGene.hpp"
#include "version3/data/ConnectionGene.hpp"
#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"

namespace Operator {

class NodeMutationParams;

Genome nodeMutation(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const NodeMutationParams& params);

class NodeMutationParams {
public:
    NodeMutationParams() = delete;
    NodeMutationParams(NodeGeneAttributes nodeAttributes);

protected:
    friend Genome nodeMutation(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const NodeMutationParams& params);
    
    const NodeGeneAttributes _nodeAttributes;  // Activation type for new hidden node
};

}