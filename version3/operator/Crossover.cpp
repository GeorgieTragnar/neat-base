#include "Crossover.hpp"
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

namespace Operator {

namespace {
    thread_local std::random_device rd;
    thread_local std::mt19937 gen(rd());
    thread_local std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    struct NodeGeneInfo {
        uint32_t innovationNumber;
        const NodeGene* genePtr;
        bool fromParentA;
    };

    struct ConnectionGeneInfo {
        uint32_t innovationNumber;
        const ConnectionGene* genePtr;
        const void* rawDataPtr;
        bool fromParentA;
    };

    struct InnovationRange {
        uint32_t min;
        uint32_t max;
    };

    InnovationRange getInnovationRange(const std::vector<ConnectionGene>& connections) {
        if (connections.empty()) {
            return {0, 0};
        }
        
        uint32_t min = connections[0].get_historyID();
        uint32_t max = connections[0].get_historyID();
        
        for (const auto& conn : connections) {
            uint32_t id = conn.get_historyID();
            if (id < min) min = id;
            if (id > max) max = id;
        }
        
        return {min, max};
    }
}

Genome crossover(
    const Genome& parentA, 
    const Analysis::GenomeAnalytics& analyticsA,
    const Genome& parentB, 
    const Analysis::GenomeAnalytics& analyticsB,
    const CrossoverParams& params) {
    
    // Input validation
    const auto& nodesA = parentA.get_nodeGenes();
    const auto& nodesB = parentB.get_nodeGenes();
    const auto& connectionsA = parentA.get_connectionGenes();
    const auto& connectionsB = parentB.get_connectionGenes();
    
    // Assert I/O structure compatibility
    assert(!nodesA.empty() && !nodesB.empty());
    
    // Determine fitter parent
    bool parentAIsFitter = analyticsA.getFitness().isBetterThan(analyticsB.getFitness());
    bool equalFitness = analyticsA.getFitness().isEqualTo(analyticsB.getFitness());
    
    if (equalFitness) {
        // For equal fitness, randomly choose which parent contributes disjoint/excess genes
        parentAIsFitter = uniform_dist(gen) < 0.5;
    }
    
    // Create maps for efficient lookup
    std::unordered_map<uint32_t, const NodeGene*> nodeMapA, nodeMapB;
    std::unordered_map<uint32_t, const ConnectionGene*> connMapA, connMapB;
    
    for (const auto& node : nodesA) {
        nodeMapA[node.get_historyID()] = &node;
    }
    for (const auto& node : nodesB) {
        nodeMapB[node.get_historyID()] = &node;
    }
    for (const auto& conn : connectionsA) {
        connMapA[conn.get_historyID()] = &conn;
    }
    for (const auto& conn : connectionsB) {
        connMapB[conn.get_historyID()] = &conn;
    }

    // Categorize node genes
    std::vector<NodeGeneInfo> matchingNodes, disjointExcessNodes;
    std::unordered_set<uint32_t> processedNodes;
    
    for (const auto& node : nodesA) {
        uint32_t id = node.get_historyID();
        processedNodes.insert(id);
        
        if (nodeMapB.find(id) != nodeMapB.end()) {
            matchingNodes.push_back({id, &node, true});
        } else {
            disjointExcessNodes.push_back({id, &node, true});
        }
    }
    
    for (const auto& node : nodesB) {
        uint32_t id = node.get_historyID();
        if (processedNodes.find(id) == processedNodes.end()) {
            disjointExcessNodes.push_back({id, &node, false});
        }
    }

    // Categorize connection genes
    std::vector<ConnectionGeneInfo> matchingConnections, disjointExcessConnections;
    std::unordered_set<uint32_t> processedConnections;
    
    for (const auto& conn : connectionsA) {
        uint32_t id = conn.get_historyID();
        processedConnections.insert(id);
        
        if (connMapB.find(id) != connMapB.end()) {
            matchingConnections.push_back({id, &conn, conn.get_rawData(), true});
        } else {
            disjointExcessConnections.push_back({id, &conn, conn.get_rawData(), true});
        }
    }
    
    for (const auto& conn : connectionsB) {
        uint32_t id = conn.get_historyID();
        if (processedConnections.find(id) == processedConnections.end()) {
            disjointExcessConnections.push_back({id, &conn, conn.get_rawData(), false});
        } else {
            // This is a matching connection from parent B
            matchingConnections.push_back({id, &conn, conn.get_rawData(), false});
        }
    }
    
    // Select inherited node genes
    std::vector<const NodeGene*> selectedNodes;
    
    // Matching nodes: always include (could randomly select between parents, but both should be identical)
    for (const auto& nodeInfo : matchingNodes) {
        selectedNodes.push_back(nodeInfo.genePtr);
    }
    
    // Disjoint/excess nodes: only from fitter parent
    for (const auto& nodeInfo : disjointExcessNodes) {
        if ((parentAIsFitter && nodeInfo.fromParentA) || 
            (!parentAIsFitter && !nodeInfo.fromParentA)) {
            selectedNodes.push_back(nodeInfo.genePtr);
        }
    }
    
    // Select inherited connection genes
    std::vector<const void*> selectedConnections;
    std::vector<uint32_t> connectionsToReactivate;
    
    // Matching connections: random selection (50/50) between parents
    std::unordered_map<uint32_t, std::vector<ConnectionGeneInfo>> matchingByID;
    for (const auto& connInfo : matchingConnections) {
        matchingByID[connInfo.innovationNumber].push_back(connInfo);
    }
    
    for (const auto& [innovationID, connections] : matchingByID) {
        // Should have exactly 2 connections (one from each parent)
        assert(connections.size() == 2);
        
        // Randomly select from parent A or B
        bool selectFromA = uniform_dist(gen) < 0.5;
        const auto& selectedConn = selectFromA ? 
            (connections[0].fromParentA ? connections[0] : connections[1]) :
            (connections[0].fromParentA ? connections[1] : connections[0]);
        
        selectedConnections.push_back(selectedConn.rawDataPtr);
        
        // Apply disabled gene reactivation probability
        if (!reinterpret_cast<const ConnectionGene::RawData*>(selectedConn.rawDataPtr)->_attributes.enabled) {
            if (uniform_dist(gen) < params._disabledGeneReactivationProbability) {
                connectionsToReactivate.push_back(selectedConn.innovationNumber);
            }
        }
    }
    
    // Disjoint/excess connections: only from fitter parent
    for (const auto& connInfo : disjointExcessConnections) {
        if ((parentAIsFitter && connInfo.fromParentA) || 
            (!parentAIsFitter && !connInfo.fromParentA)) {
            selectedConnections.push_back(connInfo.rawDataPtr);
        }
    }
    
    // Build offspring
    RawGenomeParams rawParams;
    rawParams._nodeGenes = selectedNodes;
    rawParams._rawConnectionGeneData = selectedConnections;
    
    Genome offspring(rawParams);
    
    // Apply disabled gene reactivation to the offspring
    for (uint32_t historyID : connectionsToReactivate) {
        auto& connections = offspring.get_connectionGenes();
        for (auto& conn : connections) {
            if (conn.get_historyID() == historyID) {
                conn.get_attributes().enabled = true;
                break;
            }
        }
    }
    
    return offspring;
}

}