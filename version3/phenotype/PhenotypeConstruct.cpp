#include "PhenotypeConstruct.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <climits>

namespace Operator {

void phenotypeConstruct(Genome& genome)
{
	auto& phenotype = genome.get_phenotype();
	auto& connectionDeltas = genome.get_connectionGeneDeltas();
	auto& nodeDeltas = genome.get_nodeGeneDeltas();
	
	phenotype._inputIndices.clear();
	phenotype._outputIndices.clear();
	phenotype._nodeGeneAttributes.clear();
	phenotype._orderedConnections.clear();
	phenotype._biasIndex = SIZE_MAX; // Initialize to invalid index

	std::unordered_set<uint32_t> includedNodeGeneHistoryIDs;

	for (const auto& node : genome.get_nodeGenes()) {
		if (node.get_type() == NodeType::INPUT || 
			node.get_type() == NodeType::OUTPUT || 
			node.get_type() == NodeType::BIAS) {
			includedNodeGeneHistoryIDs.insert(node.get_historyID());
		}
	}

	const auto& nodes = genome.get_nodeGenes();
	for (const auto& conn : genome.get_connectionGenes()) {
		if (conn.get_attributes().enabled) {
			includedNodeGeneHistoryIDs.insert(nodes[conn.get_sourceNodeIndex()].get_historyID());
			includedNodeGeneHistoryIDs.insert(nodes[conn.get_targetNodeIndex()].get_historyID());
		}
	}

	std::unordered_map<uint32_t, size_t> historyIDToIndex;
	phenotype._nodeGeneAttributes.reserve(includedNodeGeneHistoryIDs.size());

	size_t nodeIndex = 0;
	for (const auto& node : genome.get_nodeGenes()) {
		if (includedNodeGeneHistoryIDs.count(node.get_historyID()) > 0) {
			historyIDToIndex[node.get_historyID()] = nodeIndex;
			phenotype._nodeGeneAttributes.push_back(node.get_attributes());

			if (node.get_type() == NodeType::INPUT) {
				phenotype._inputIndices.push_back(nodeIndex);
			} else if (node.get_type() == NodeType::OUTPUT) {
				phenotype._outputIndices.push_back(nodeIndex);
			} else if (node.get_type() == NodeType::BIAS) {
				phenotype._biasIndex = nodeIndex;
			}
			
			nodeIndex++;
		}
	}

	size_t enabledConnectionCount = 0;
	for (const auto& conn : genome.get_connectionGenes()) {
		if (conn.get_attributes().enabled) {
			enabledConnectionCount++;
		}
	}

	phenotype._orderedConnections.reserve(enabledConnectionCount);

	for (const auto& conn : genome.get_connectionGenes()) {
		if (conn.get_attributes().enabled) {
			uint32_t sourceID = nodes[conn.get_sourceNodeIndex()].get_historyID();
			uint32_t targetID = nodes[conn.get_targetNodeIndex()].get_historyID();

			Phenotype::Connection phenConn;
			phenConn._sourceNodeIndex = historyIDToIndex[sourceID];
			phenConn._targetNodeIndex = historyIDToIndex[targetID];
			phenConn._connectionGeneAttribute = conn.get_attributes();

			phenotype._orderedConnections.push_back(std::move(phenConn));
		}
	}
	
	// Clear both delta vectors after reconstruction
	connectionDeltas.clear();
	nodeDeltas.clear();
}

}