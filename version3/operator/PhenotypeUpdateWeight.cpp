#include "PhenotypeUpdateWeight.hpp"
#include <unordered_map>
#include <cassert>

namespace Operator {

void phenotypeUpdateWeight(Genome& genome)
{
	auto& phenotype = genome.get_phenotype();
	auto& connectionDeltas = genome.get_connectionGeneDeltas();
	const auto& connectionGenes = genome.get_connectionGenes();
	
	// Assert that connection deltas is non-empty
	assert(!connectionDeltas.empty() && "Connection deltas must be non-empty for weight phenotype update");
	
	// Create a map of history ID to connection gene for quick lookup
	std::unordered_map<uint32_t, const ConnectionGene*> historyIDToConnection;
	for (const auto& conn : connectionGenes) {
		historyIDToConnection[conn.get_historyID()] = &conn;
	}
	
	// Create a map of history ID pairs to phenotype connection index for quick lookup
	std::unordered_map<uint64_t, size_t> connectionPairToIndex;
	for (size_t i = 0; i < phenotype._orderedConnections.size(); ++i) {
		const auto& phenConn = phenotype._orderedConnections[i];
		
		// Find the corresponding genome nodes to get their history IDs
		size_t sourceNodeIndex = phenConn._sourceNodeIndex;
		size_t targetNodeIndex = phenConn._targetNodeIndex;
		
		// We need to map from phenotype node indices back to history IDs
		// This requires iterating through genome nodes to find matches
		uint32_t sourceHistoryID = 0;
		uint32_t targetHistoryID = 0;
		
		size_t nodeIndex = 0;
		for (const auto& node : genome.get_nodeGenes()) {
			if (nodeIndex == sourceNodeIndex) {
				sourceHistoryID = node.get_historyID();
			}
			if (nodeIndex == targetNodeIndex) {
				targetHistoryID = node.get_historyID();
			}
			nodeIndex++;
			if (sourceHistoryID != 0 && targetHistoryID != 0) break;
		}
		
		// Create a unique key from source and target history IDs
		uint64_t connectionKey = (static_cast<uint64_t>(sourceHistoryID) << 32) | targetHistoryID;
		connectionPairToIndex[connectionKey] = i;
	}
	
	// Update phenotype weights using genome as source of truth
	for (uint32_t historyID : connectionDeltas) {
		auto connectionIt = historyIDToConnection.find(historyID);
		assert(connectionIt != historyIDToConnection.end() && "Connection with history ID not found in genome");
		
		const ConnectionGene* connection = connectionIt->second;
		
		// Only update if connection is enabled
		if (connection->get_attributes().enabled) {
			uint32_t sourceHistoryID = connection->get_sourceNodeGene().get_historyID();
			uint32_t targetHistoryID = connection->get_targetNodeGene().get_historyID();
			
			uint64_t connectionKey = (static_cast<uint64_t>(sourceHistoryID) << 32) | targetHistoryID;
			auto indexIt = connectionPairToIndex.find(connectionKey);
			
			if (indexIt != connectionPairToIndex.end()) {
				// Update the weight in the phenotype
				size_t phenotypeIndex = indexIt->second;
				phenotype._orderedConnections[phenotypeIndex]._connectionGeneAttribute.weight = 
					connection->get_attributes().weight;
			}
		}
	}
	
	// Clear connection deltas after update
	connectionDeltas.clear();
}

}