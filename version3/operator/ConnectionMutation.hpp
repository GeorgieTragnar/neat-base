#pragma once

#include <vector>
#include <memory>

#include "version3/data/NodeGene.hpp"
#include "version3/data/ConnectionGene.hpp"
#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"

namespace Operator {

class ConnectionMutationParams;

Genome connectionMutation(const Genome& genome, std::unique_ptr<HistoryTracker> historyTracker, const ConnectionMutationParams& params);

class ConnectionMutationParams {
public:
    enum class NetworkTopology {
        FEED_FORWARD,    // Topology hint for analysis stage
        RECURRENT        // Topology hint for analysis stage
    };
    
    ConnectionMutationParams() = delete;
    ConnectionMutationParams(double connectionRate,
                           double weightRange,
                           NetworkTopology topology);

protected:
    friend Genome connectionMutation(const Genome& genome, std::unique_ptr<HistoryTracker> historyTracker, const ConnectionMutationParams& params);
    
    const double _connectionRate;        // Probability to add a connection [0.0, 1.0]
    const double _weightRange;           // Range for new connection weights [-range, +range] (> 0)
    const NetworkTopology _topology;    // Topology hint for downstream analysis
};

}