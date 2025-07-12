#include "HasPossibleConnections.hpp"

namespace Operator {

namespace {
    // Helper to check if a connection already exists between two nodes
    bool connectionExists(const std::vector<ConnectionGene>& connections,
                         const std::vector<NodeGene>& nodes,
                         uint32_t sourceHistoryID, 
                         uint32_t targetHistoryID) {
        for (const auto& conn : connections) {
            uint32_t connSourceID = nodes[conn.get_sourceNodeIndex()].get_historyID();
            uint32_t connTargetID = nodes[conn.get_targetNodeIndex()].get_historyID();
            if (connSourceID == sourceHistoryID && connTargetID == targetHistoryID) {
                return true;
            }
        }
        return false;
    }
}

bool hasPossibleConnections(const Genome& genome) {
    // Access genome data through friend access
    const auto& nodes = const_cast<Genome&>(genome).get_nodeGenes();
    const auto& connections = const_cast<Genome&>(genome).get_connectionGenes();
    
    // Need at least 2 nodes to create a connection
    if (nodes.size() < 2) {
        return false;
    }
    
    // Check all possible node pairs for valid connections
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
            if (connectionExists(connections, nodes, sourceNode.get_historyID(), targetNode.get_historyID())) {
                continue;
            }
            
            // Found a possible connection
            return true;
        }
    }
    
    // No possible connections found
    return false;
}

} // namespace Operator