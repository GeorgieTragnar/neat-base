#include "PhenotypeUpdateConnection.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cassert>

namespace Operator {

void phenotypeUpdateConnection(Genome& genome)
{
	auto& phenotype = genome.get_phenotype();
	auto& connectionDeltas = genome.get_connectionGeneDeltas();
	const auto& connectionGenes = genome.get_connectionGenes();
	const auto& nodeGenes = genome.get_nodeGenes();
	
	// Early return if no connection deltas to process (valid state - no mutations occurred)
	if (connectionDeltas.empty()) {
		return;
	}
	
	// Assert that connection deltas contains exactly 1 ID when mutations did occur
	assert(connectionDeltas.size() == 1 && "Connection deltas must contain exactly 1 ID for connection phenotype update");
	
	uint32_t newConnectionHistoryID = connectionDeltas[0];
	
	// Find the new connection in the genome
	const ConnectionGene* newConnection = nullptr;
	for (const auto& conn : connectionGenes) {
		if (conn.get_historyID() == newConnectionHistoryID) {
			newConnection = &conn;
			break;
		}
	}
	assert(newConnection != nullptr && "New connection with history ID not found in genome");
	assert(newConnection->get_attributes().enabled && "New connection must be enabled");
	
	uint32_t sourceHistoryID = nodeGenes[newConnection->get_sourceNodeIndex()].get_historyID();
	uint32_t targetHistoryID = nodeGenes[newConnection->get_targetNodeIndex()].get_historyID();
	
	// Build current phenotype node mapping (history ID to phenotype index)
	std::unordered_map<uint32_t, size_t> historyIDToIndex;
	std::unordered_set<uint32_t> includedNodeHistoryIDs;
	
	// First, determine which nodes are currently included by scanning the ENTIRE genome state
	// (as if we were doing a full reconstruction, but excluding the new connection)
	for (const auto& node : nodeGenes) {
		if (node.get_type() == NodeType::INPUT || 
			node.get_type() == NodeType::OUTPUT || 
			node.get_type() == NodeType::BIAS) {
			includedNodeHistoryIDs.insert(node.get_historyID());
		}
	}
	
	for (const auto& conn : connectionGenes) {
		if (conn.get_attributes().enabled && conn.get_historyID() != newConnectionHistoryID) {
			includedNodeHistoryIDs.insert(nodeGenes[conn.get_sourceNodeIndex()].get_historyID());
			includedNodeHistoryIDs.insert(nodeGenes[conn.get_targetNodeIndex()].get_historyID());
		}
	}
	
	// Check if we need to add new nodes to phenotype
	bool needToAddSourceNode = includedNodeHistoryIDs.find(sourceHistoryID) == includedNodeHistoryIDs.end();
	bool needToAddTargetNode = includedNodeHistoryIDs.find(targetHistoryID) == includedNodeHistoryIDs.end();
	
	if (needToAddSourceNode || needToAddTargetNode) {
		// Need to rebuild phenotype node list to include new nodes
		phenotype._nodeGeneAttributes.clear();
		phenotype._inputIndices.clear();
		phenotype._outputIndices.clear();
		phenotype._orderedConnections.clear(); // Clear existing connections too!
		
		// Add the new connection nodes to included set
		includedNodeHistoryIDs.insert(sourceHistoryID);
		includedNodeHistoryIDs.insert(targetHistoryID);
		
		size_t nodeIndex = 0;
		for (const auto& node : nodeGenes) {
			if (includedNodeHistoryIDs.count(node.get_historyID()) > 0) {
				historyIDToIndex[node.get_historyID()] = nodeIndex;
				phenotype._nodeGeneAttributes.push_back(node.get_attributes());
				
				if (node.get_type() == NodeType::INPUT) {
					phenotype._inputIndices.push_back(nodeIndex);
				} else if (node.get_type() == NodeType::OUTPUT) {
					phenotype._outputIndices.push_back(nodeIndex);
				}
				
				nodeIndex++;
			}
		}
		
		// Rebuild all enabled connections with new indices
		for (const auto& conn : connectionGenes) {
			if (conn.get_attributes().enabled) {
				uint32_t connSourceID = nodeGenes[conn.get_sourceNodeIndex()].get_historyID();
				uint32_t connTargetID = nodeGenes[conn.get_targetNodeIndex()].get_historyID();
				
				Phenotype::Connection phenConn;
				phenConn._sourceNodeIndex = historyIDToIndex[connSourceID];
				phenConn._targetNodeIndex = historyIDToIndex[connTargetID];
				phenConn._connectionGeneAttribute = conn.get_attributes();
				
				phenotype._orderedConnections.push_back(std::move(phenConn));
			}
		}
	} else {
		// Nodes already exist, just build the mapping
		size_t nodeIndex = 0;
		for (const auto& node : nodeGenes) {
			if (includedNodeHistoryIDs.count(node.get_historyID()) > 0) {
				historyIDToIndex[node.get_historyID()] = nodeIndex;
				nodeIndex++;
			}
		}
		
		// Add the new connection to phenotype
		Phenotype::Connection phenConn;
		phenConn._sourceNodeIndex = historyIDToIndex[sourceHistoryID];
		phenConn._targetNodeIndex = historyIDToIndex[targetHistoryID];
		phenConn._connectionGeneAttribute = newConnection->get_attributes();
		
		phenotype._orderedConnections.push_back(std::move(phenConn));
	}
	
	// Clear connection deltas after update
	connectionDeltas.clear();
}

}