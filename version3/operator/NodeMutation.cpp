#include "NodeMutation.hpp"
#include <random>
#include <cassert>
#include <algorithm>

using namespace Operator;

NodeMutationParams::NodeMutationParams(NodeGeneAttributes nodeAttributes)
    : _nodeAttributes(nodeAttributes)
{
    // No parameter validation needed - NodeGeneAttributes is always valid
}

namespace {
    // Helper to find all enabled connection indices
    std::vector<size_t> findEnabledConnectionIndices(const std::vector<ConnectionGene>& connections) {
        std::vector<size_t> enabledIndices;
        
        for (size_t i = 0; i < connections.size(); ++i) {
            if (connections[i].get_attributes().enabled) {
                enabledIndices.push_back(i);
            }
        }
        
        return enabledIndices;
    }
    
    // Helper to check if a node with given history ID exists in the genome
    bool nodeExistsInGenome(const std::vector<NodeGene>& nodes, uint32_t historyID) {
        for (const auto& node : nodes) {
            if (node.get_historyID() == historyID) {
                return true;
            }
        }
        return false;
    }
    
    // Helper to check if a node has active (enabled) connections in the genome
    bool nodeHasActiveConnections(const std::vector<ConnectionGene>& connections, 
                                 const std::vector<NodeGene>& nodes, 
                                 uint32_t nodeHistoryID) {
        for (const auto& conn : connections) {
            if (conn.get_attributes().enabled) {
                uint32_t sourceHistoryID = nodes[conn.get_sourceNodeIndex()].get_historyID();
                uint32_t targetHistoryID = nodes[conn.get_targetNodeIndex()].get_historyID();
                if (sourceHistoryID == nodeHistoryID || targetHistoryID == nodeHistoryID) {
                    return true;
                }
            }
        }
        return false;
    }
    
    // Helper to find first available node ID from a list
    uint32_t findFirstAvailableNode(const std::vector<uint32_t>& nodeIDs, 
                                   const std::vector<NodeGene>& genomeNodes) {
        for (uint32_t nodeID : nodeIDs) {
            if (!nodeExistsInGenome(genomeNodes, nodeID)) {
                return nodeID;
            }
        }
        return 0; // No available node found
    }
    
    // Helper to determine the split node ID using the innovation decision tree
    uint32_t determineSplitNodeID(uint32_t connectionID,
                                 const std::vector<NodeGene>& genomeNodes,
                                 const std::vector<ConnectionGene>& genomeConnections,
                                 std::shared_ptr<HistoryTracker> historyTracker,
                                 uint32_t sourceNodeID,
                                 uint32_t targetNodeID) {
        
        // Step 1: Get primary split node ID for this connection
        uint32_t primarySplitNodeID = historyTracker->get_splitNode(connectionID);
        
        // Step 2: Check if split node exists in genome
        if (!nodeExistsInGenome(genomeNodes, primarySplitNodeID)) {
            // Node doesn't exist - use primary node ID
            return primarySplitNodeID;
        }
        
        // Step 3: Check if primary split node is valid (not source or target of connection being split)
        if (primarySplitNodeID == sourceNodeID || primarySplitNodeID == targetNodeID) {
            // Primary node conflicts with connection being split - create new split branch immediately
            return historyTracker->create_splitBranch(connectionID);
        }
        
        // Step 4: Primary node is valid - check if it has active connections
        if (!nodeHasActiveConnections(genomeConnections, genomeNodes, primarySplitNodeID)) {
            // Node exists but has no active connections - reuse it
            return primarySplitNodeID;
        }
        
        // Step 5: Node exists with active connections - get all split node alternatives
        std::vector<uint32_t> allSplitNodes = historyTracker->get_allSplitNodes(connectionID);
        
        // Step 6: Find first available alternative that doesn't conflict
        for (uint32_t nodeID : allSplitNodes) {
            if (nodeID != sourceNodeID && nodeID != targetNodeID && 
                !nodeExistsInGenome(genomeNodes, nodeID)) {
                return nodeID;
            }
        }
        
        // Step 7: No alternatives available - create new split branch
        return historyTracker->create_splitBranch(connectionID);
    }
}

Genome Operator::nodeMutation(const Genome& genome, 
                             std::shared_ptr<HistoryTracker> historyTracker, 
                             const NodeMutationParams& params) {
    // Create a copy of the genome
    Genome mutatedGenome = genome;
    
    // Assert that deltas are empty before operation
    assert(mutatedGenome.get_connectionGeneDeltas().empty() && "Connection deltas must be empty before node mutation");
    assert(mutatedGenome.get_nodeGeneDeltas().empty() && "Node deltas must be empty before node mutation");
    
    auto& nodes = mutatedGenome.get_nodeGenes();
    auto& connections = mutatedGenome.get_connectionGenes();
    auto& connectionDeltas = mutatedGenome.get_connectionGeneDeltas();
    
    // Find all enabled connections
    auto enabledIndices = findEnabledConnectionIndices(connections);
    
    // Assert if no enabled connections available (user error - shouldn't call operator in this state)
    assert(!enabledIndices.empty());
    
    // Setup random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> connDist(0, enabledIndices.size() - 1);
    
    // Select random enabled connection to split
    size_t selectedIndex = enabledIndices[connDist(gen)];
    
    // Ensure capacity to prevent reallocations: 1 potential new node + 2 new connections
    mutatedGenome.ensureCapacity(1, 2);
    
    // Now safe to get references since no reallocations will occur
    const ConnectionGene& connectionToSplit = connections[selectedIndex];
    uint32_t connectionToSplitID = connectionToSplit.get_historyID();
    uint32_t sourceNodeID = nodes[connectionToSplit.get_sourceNodeIndex()].get_historyID();
    uint32_t targetNodeID = nodes[connectionToSplit.get_targetNodeIndex()].get_historyID();
    float originalWeight = connectionToSplit.get_attributes().weight;
    
    // Determine which split node to use via innovation decision tree
    uint32_t splitNodeID = determineSplitNodeID(
        connectionToSplitID,
        nodes,
        connections,
        historyTracker,
        sourceNodeID,
        targetNodeID
    );
    
    // Check if we need to create a new node or reuse existing one
    bool createNewNode = !nodeExistsInGenome(nodes, splitNodeID);
    
    if (createNewNode) {
        // Create new hidden node (safe - capacity ensured)
        nodes.emplace_back(splitNodeID, NodeType::HIDDEN, params._nodeAttributes);
    }
    
    // Find the split node in the genome (either newly created or existing)
    const NodeGene* splitNode = nullptr;
    for (const auto& node : nodes) {
        if (node.get_historyID() == splitNodeID) {
            splitNode = &node;
            break;
        }
    }
    assert(splitNode != nullptr); // Should always find the node
    
    // Get innovation numbers for new connections
    uint32_t inputConnID = historyTracker->get_connection(sourceNodeID, splitNodeID);
    uint32_t outputConnID = historyTracker->get_connection(splitNodeID, targetNodeID);
    
    // Find the indices of source, target, and split nodes in the genome
    size_t sourceNodeIndex = SIZE_MAX;
    size_t targetNodeIndex = SIZE_MAX;
    size_t splitNodeIndex = SIZE_MAX;
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        uint32_t nodeID = nodes[i].get_historyID();
        if (nodeID == sourceNodeID) {
            sourceNodeIndex = i;
        } else if (nodeID == targetNodeID) {
            targetNodeIndex = i;
        } else if (nodeID == splitNodeID) {
            splitNodeIndex = i;
        }
    }
    assert(sourceNodeIndex != SIZE_MAX && targetNodeIndex != SIZE_MAX && splitNodeIndex != SIZE_MAX);
    
    // Disable the original connection (safe - no reallocation)
    ConnectionGeneAttributes& originalConnAttribs = connections[selectedIndex].get_attributes();
    originalConnAttribs.enabled = false;
    
    // Add new connections (safe - capacity ensured)
    // Input connection: source → splitNode (weight = 1.0)
    ConnectionGeneAttributes inputConnAttribs;
    inputConnAttribs.weight = 1.0f;
    inputConnAttribs.enabled = true;
    connections.emplace_back(inputConnID, sourceNodeIndex, splitNodeIndex, inputConnAttribs);
    
    // Output connection: splitNode → target (weight = original)
    ConnectionGeneAttributes outputConnAttribs;
    outputConnAttribs.weight = originalWeight;
    outputConnAttribs.enabled = true;
    connections.emplace_back(outputConnID, splitNodeIndex, targetNodeIndex, outputConnAttribs);
    
    // Populate connection deltas: 1 disabled + 2 added (exactly 3 as per TODO spec)
    connectionDeltas.push_back(connectionToSplitID);  // disabled connection
    connectionDeltas.push_back(inputConnID);          // first added connection
    connectionDeltas.push_back(outputConnID);         // second added connection
    
    return mutatedGenome;
}