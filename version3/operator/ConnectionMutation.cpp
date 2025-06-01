#include "ConnectionMutation.hpp"
#include <random>
#include <cassert>
#include <algorithm>

using namespace Operator;

ConnectionMutationParams::ConnectionMutationParams(double connectionRate,
                                                 double weightRange,
                                                 NetworkTopology topology)
    : _connectionRate(connectionRate)
    , _weightRange(weightRange)
    , _topology(topology)
{
    // Parameter validation
    assert(connectionRate > 0.0 && connectionRate <= 1.0);
    assert(weightRange > 0.0);
}

namespace {
    // Helper to check if a connection already exists between two nodes
    bool connectionExists(const std::vector<ConnectionGene>& connections, 
                         uint32_t sourceHistoryID, 
                         uint32_t targetHistoryID) {
        for (const auto& conn : connections) {
            if (conn.get_sourceNodeGene().get_historyID() == sourceHistoryID &&
                conn.get_targetNodeGene().get_historyID() == targetHistoryID) {
                return true;
            }
        }
        return false;
    }
    
    // Helper to find all possible connection pairs (existing connections excluded)
    std::vector<std::pair<size_t, size_t>> findPossibleConnectionPairs(
        const std::vector<NodeGene>& nodes,
        const std::vector<ConnectionGene>& connections) {
        
        std::vector<std::pair<size_t, size_t>> possiblePairs;
        
        for (size_t sourceIdx = 0; sourceIdx < nodes.size(); ++sourceIdx) {
            for (size_t targetIdx = 0; targetIdx < nodes.size(); ++targetIdx) {
                const NodeGene& sourceNode = nodes[sourceIdx];
                const NodeGene& targetNode = nodes[targetIdx];
                
                // Skip self-connections
                if (sourceIdx == targetIdx) {
                    continue;
                }
                
                // Mechanical connection rules: OUTPUT nodes cannot be sources
                if (sourceNode.get_type() == NodeType::OUTPUT) {
                    continue;
                }
                
                // Mechanical connection rules: INPUT nodes cannot be targets
                if (targetNode.get_type() == NodeType::INPUT) {
                    continue;
                }
                
                // Check if connection already exists
                if (connectionExists(connections, sourceNode.get_historyID(), targetNode.get_historyID())) {
                    continue;
                }
                
                possiblePairs.push_back({sourceIdx, targetIdx});
            }
        }
        
        return possiblePairs;
    }
}

Genome Operator::connectionMutation(const Genome& genome, 
                                   std::unique_ptr<HistoryTracker> historyTracker, 
                                   const ConnectionMutationParams& params) {
    // Create a copy of the genome
    Genome mutatedGenome = genome;
    
    const auto& nodes = mutatedGenome.get_nodeGenes();
    auto& connections = mutatedGenome.get_connectionGenes();
    
    // Assert if insufficient nodes to create connection
    assert(nodes.size() >= 2);
    
    // Find all possible connection pairs
    auto possiblePairs = findPossibleConnectionPairs(nodes, connections);
    
    // Assert if no connections possible (user error - shouldn't call operator in this state)
    assert(!possiblePairs.empty());
    
    // Setup random number generators
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> weightDist(-params._weightRange, params._weightRange);
    
    // Select a random possible pair
    std::uniform_int_distribution<size_t> pairDist(0, possiblePairs.size() - 1);
    auto selectedPair = possiblePairs[pairDist(gen)];
    
    const NodeGene& sourceNode = nodes[selectedPair.first];
    const NodeGene& targetNode = nodes[selectedPair.second];
    
    // Get innovation number for this connection
    uint32_t innovationNumber = historyTracker->get_connection(
        sourceNode.get_historyID(), 
        targetNode.get_historyID()
    );
    
    // Create new connection attributes
    ConnectionGeneAttributes newConnAttribs;
    newConnAttribs.weight = static_cast<float>(weightDist(gen));
    newConnAttribs.enabled = true;
    
    // Create new connections vector to completely avoid reallocation issues
    std::vector<ConnectionGene> newConnections;
    newConnections.reserve(connections.size() + 1);
    
    // Copy all existing connections
    for (const auto& connection : connections) {
        ConnectionGeneAttributes attrs = connection.get_attributes();
        newConnections.emplace_back(connection.get_historyID(),
                                   connection.get_sourceNodeGene(),
                                   connection.get_targetNodeGene(),
                                   attrs);
    }
    
    // Add the new connection
    newConnections.emplace_back(innovationNumber, sourceNode, targetNode, newConnAttribs);
    
    // Replace the original connections vector
    connections = std::move(newConnections);
    
    return mutatedGenome;
}