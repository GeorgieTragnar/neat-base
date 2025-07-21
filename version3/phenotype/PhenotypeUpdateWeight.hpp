#pragma once

#include "version3/data/Genome.hpp"
#include "version3/data/PopulationContainer.hpp"
#include <cstdint>
#include <cassert>

namespace Operator {

// Container-based operator
template<typename FitnessResultType>
void phenotypeUpdateWeight(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& generation
);

// Template implementation
template<typename FitnessResultType>
void phenotypeUpdateWeight(
    PopulationContainer<FitnessResultType>& container,
    const size_t& targetIndex,
    const uint32_t& generation
) {
    // Get mutable access to genomes at specified generation
    auto& genomes = container.getCurrentGenomes(generation);
    
    // Validate bounds
    assert(targetIndex < genomes.size() && 
           "phenotypeUpdateWeight: targetIndex out of bounds");
    
    auto& genome = genomes[targetIndex];
	auto& phenotype = genome.get_phenotype();
	auto& connectionDeltas = genome.get_connectionGeneDeltas();
	const auto& connectionGenes = genome.get_connectionGenes();
	
	// static auto logger = LOGGER("operator.PhenotypeUpdateWeight");
	// LOG_DEBUG("phenotypeUpdateWeight called with {} connection deltas", connectionDeltas.size());
	
	// Early return if no connection deltas to process (valid state - no mutations occurred)
	if (connectionDeltas.empty()) {
		// LOG_DEBUG("phenotypeUpdateWeight: no deltas to process, returning early");
		return;
	}
	
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
			const auto& nodeGenes = genome.get_nodeGenes();
			uint32_t sourceHistoryID = nodeGenes[connection->get_sourceNodeIndex()].get_historyID();
			uint32_t targetHistoryID = nodeGenes[connection->get_targetNodeIndex()].get_historyID();
			
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