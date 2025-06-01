#pragma once

#include <vector>
#include <memory>

#include "../data/NodeGene.hpp"
#include "../data/ConnectionGene.hpp"
#include "../data/Genome.hpp"
#include "../data/HistoryTracker.hpp"

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