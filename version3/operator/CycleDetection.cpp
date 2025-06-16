#include "CycleDetection.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cassert>

using namespace Operator;

// Parameter validation constructor
CycleDetectionParams::CycleDetectionParams(bool includeAllNodeTypes, bool failFast)
    : _includeAllNodeTypes(includeAllNodeTypes), _failFast(failFast) {
    // No parameter validation needed - all boolean values are valid
}

namespace {
    // Node visit states for DFS cycle detection
    enum class VisitState {
        UNVISITED,  // Node not yet visited
        VISITING,   // Node currently in DFS path (grey)
        VISITED     // Node completely processed (black)
    };
    
    // DFS-based cycle detection with path tracking
    bool dfsHasCycle(const std::unordered_map<uint32_t, std::vector<uint32_t>>& adjacencyList,
                     uint32_t currentNode, 
                     std::unordered_map<uint32_t, VisitState>& visitStates,
                     bool failFast) {
        
        // Mark current node as being visited (in current path)
        visitStates[currentNode] = VisitState::VISITING;
        
        // Check all neighbors
        auto it = adjacencyList.find(currentNode);
        if (it != adjacencyList.end()) {
            for (uint32_t neighbor : it->second) {
                VisitState neighborState = visitStates[neighbor];
                
                if (neighborState == VisitState::VISITING) {
                    // Back edge found - cycle detected!
                    return true;
                } else if (neighborState == VisitState::UNVISITED) {
                    // Continue DFS
                    if (dfsHasCycle(adjacencyList, neighbor, visitStates, failFast)) {
                        return true;
                    }
                }
                // If neighbor is VISITED, it's already been completely processed - safe to ignore
            }
        }
        
        // Mark current node as completely processed
        visitStates[currentNode] = VisitState::VISITED;
        return false;
    }
}

bool Operator::hasCycles(const Genome& genome, const CycleDetectionParams& params) {
    const auto& nodes = genome.get_nodeGenes();
    const auto& connections = genome.get_connectionGenes();
    
    // Build graph representation
    std::unordered_map<uint32_t, std::vector<uint32_t>> adjacencyList;
    std::unordered_set<uint32_t> graphNodes;
    
    // Add all relevant nodes to graph (even if they have no connections)
    for (const auto& node : nodes) {
        uint32_t nodeId = node.get_historyID();
        NodeType nodeType = node.get_type();
        
        // Include node based on configuration
        if (params._includeAllNodeTypes || nodeType == NodeType::HIDDEN) {
            graphNodes.insert(nodeId);
            adjacencyList[nodeId] = std::vector<uint32_t>();
        }
    }
    
    // Add edges from enabled connections only
    for (const auto& conn : connections) {
        const ConnectionGeneAttributes& attrs = conn.get_attributes();
        
        // Skip disabled connections - they don't affect execution flow
        if (!attrs.enabled) {
            continue;
        }
        
        uint32_t sourceId = nodes[conn.get_sourceNodeIndex()].get_historyID();
        uint32_t targetId = nodes[conn.get_targetNodeIndex()].get_historyID();
        
        // Only add edge if both nodes are included in our analysis
        if (graphNodes.find(sourceId) != graphNodes.end() && 
            graphNodes.find(targetId) != graphNodes.end()) {
            adjacencyList[sourceId].push_back(targetId);
        }
    }
    
    // Early exit for empty graph
    if (graphNodes.empty()) {
        return false;
    }
    
    // Initialize visit states
    std::unordered_map<uint32_t, VisitState> visitStates;
    for (uint32_t nodeId : graphNodes) {
        visitStates[nodeId] = VisitState::UNVISITED;
    }
    
    // Check each disconnected component for cycles
    for (uint32_t nodeId : graphNodes) {
        if (visitStates[nodeId] == VisitState::UNVISITED) {
            if (dfsHasCycle(adjacencyList, nodeId, visitStates, params._failFast)) {
                return true; // Cycle found in this component
            }
        }
    }
    
    return false; // No cycles found in any component
}