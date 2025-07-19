#include "PhenotypeUpdateNode.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cassert>

namespace Operator {

void phenotypeUpdateNode(Genome& genome)
{
	auto& phenotype = genome.get_phenotype();
	auto& connectionDeltas = genome.get_connectionGeneDeltas();
	const auto& connectionGenes = genome.get_connectionGenes();
	const auto& nodeGenes = genome.get_nodeGenes();
	
	// Early return if no connection deltas to process (valid state - no mutations occurred)
	if (connectionDeltas.empty()) {
		return;
	}
	
	// Assert that connection deltas contains exactly 3 IDs when mutations did occur
	assert(connectionDeltas.size() == 3 && "Connection deltas must contain exactly 3 IDs for node phenotype update");
	
	uint32_t disabledConnectionID = connectionDeltas[0];  // First ID is the disabled connection
	uint32_t firstNewConnectionID = connectionDeltas[1];  // Second ID is first added connection
	uint32_t secondNewConnectionID = connectionDeltas[2]; // Third ID is second added connection
	
	// Find the connections in the genome
	const ConnectionGene* disabledConnection = nullptr;
	const ConnectionGene* firstNewConnection = nullptr;
	const ConnectionGene* secondNewConnection = nullptr;
	
	for (const auto& conn : connectionGenes) {
		uint32_t historyID = conn.get_historyID();
		if (historyID == disabledConnectionID) {
			disabledConnection = &conn;
		} else if (historyID == firstNewConnectionID) {
			firstNewConnection = &conn;
		} else if (historyID == secondNewConnectionID) {
			secondNewConnection = &conn;
		}
	}
	
	assert(disabledConnection != nullptr && "Disabled connection not found in genome");
	assert(firstNewConnection != nullptr && "First new connection not found in genome");
	assert(secondNewConnection != nullptr && "Second new connection not found in genome");
	assert(!disabledConnection->get_attributes().enabled && "Disabled connection must be disabled");
	assert(firstNewConnection->get_attributes().enabled && "First new connection must be enabled");
	assert(secondNewConnection->get_attributes().enabled && "Second new connection must be enabled");
	
	// Completely rebuild the phenotype to handle the node mutation changes
	// This is simpler and more reliable than trying to incrementally update
	
	phenotype._inputIndices.clear();
	phenotype._outputIndices.clear();
	phenotype._nodeGeneAttributes.clear();
	phenotype._orderedConnections.clear();

	std::unordered_set<uint32_t> includedNodeGeneHistoryIDs;

	// Determine which nodes should be included in the phenotype
	for (const auto& node : nodeGenes) {
		if (node.get_type() == NodeType::INPUT || 
			node.get_type() == NodeType::OUTPUT || 
			node.get_type() == NodeType::BIAS) {
			includedNodeGeneHistoryIDs.insert(node.get_historyID());
		}
	}

	for (const auto& conn : connectionGenes) {
		if (conn.get_attributes().enabled) {
			includedNodeGeneHistoryIDs.insert(nodeGenes[conn.get_sourceNodeIndex()].get_historyID());
			includedNodeGeneHistoryIDs.insert(nodeGenes[conn.get_targetNodeIndex()].get_historyID());
		}
	}

	// Build node mapping and populate phenotype nodes
	std::unordered_map<uint32_t, size_t> historyIDToIndex;
	phenotype._nodeGeneAttributes.reserve(includedNodeGeneHistoryIDs.size());

	size_t nodeIndex = 0;
	for (const auto& node : nodeGenes) {
		if (includedNodeGeneHistoryIDs.count(node.get_historyID()) > 0) {
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

	// Count enabled connections and populate phenotype connections
	size_t enabledConnectionCount = 0;
	for (const auto& conn : connectionGenes) {
		if (conn.get_attributes().enabled) {
			enabledConnectionCount++;
		}
	}

	phenotype._orderedConnections.reserve(enabledConnectionCount);

	for (const auto& conn : connectionGenes) {
		if (conn.get_attributes().enabled) {
			uint32_t sourceID = nodeGenes[conn.get_sourceNodeIndex()].get_historyID();
			uint32_t targetID = nodeGenes[conn.get_targetNodeIndex()].get_historyID();

			Phenotype::Connection phenConn;
			phenConn._sourceNodeIndex = historyIDToIndex[sourceID];
			phenConn._targetNodeIndex = historyIDToIndex[targetID];
			phenConn._connectionGeneAttribute = conn.get_attributes();

			phenotype._orderedConnections.push_back(std::move(phenConn));
		}
	}
	
	// Clear connection deltas after update
	connectionDeltas.clear();
}

}