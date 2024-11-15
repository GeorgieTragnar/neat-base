
#include "NetworkValidator.hpp"

#include "core/Genome.hpp"
#include <queue>

namespace neat {
namespace network {

// Validate all nodes exist and have valid IDs
void NetworkValidator::validateNodeStructure(const core::Genome& genome) {
	const auto& nodes = genome.getNodes();
	
	// Check node ID bounds
	for (const auto& [id, type] : nodes) {
		if (id > genome.getMaxNodeId() && id != -1) {
			throw std::runtime_error("Node ID " + std::to_string(id) + 
				" exceeds maxNodeId " + std::to_string(genome.getMaxNodeId()));
		}
	}

	// Verify we have at least one input and output node
	bool hasInput = false, hasOutput = false;
	for (const auto& [id, type] : nodes) {
		if (type == core::ENodeType::INPUT) hasInput = true;
		if (type == core::ENodeType::OUTPUT) hasOutput = true;
	}
	
	if (!hasInput || !hasOutput) {
		throw std::runtime_error("Network must have at least one input and output node");
	}
}

// Validate all genes reference valid nodes and have valid structure
void NetworkValidator::validateGeneStructure(const core::Genome& genome) {
	const auto& genes = genome.getGenes();
	const auto& nodes = genome.getNodes();
	
	for (const auto& gene : genes) {
		if (gene.enabled) {
			// Validate input node
			if (gene.inputNode != -1 && nodes.find(gene.inputNode) == nodes.end()) {
				throw std::runtime_error("Gene references non-existent input node: " + 
					std::to_string(gene.inputNode));
			}
			
			// Validate output node
			if (nodes.find(gene.outputNode) == nodes.end()) {
				throw std::runtime_error("Gene references non-existent output node: " + 
					std::to_string(gene.outputNode));
			}
			
			// Validate node types
			if (gene.inputNode != -1) {  // Skip bias node check
				auto inputType = nodes.at(gene.inputNode);
				if (inputType == core::ENodeType::OUTPUT) {
					throw std::runtime_error("Output node cannot be input to connection");
				}
			}
			
			auto outputType = nodes.at(gene.outputNode);
			if (outputType == core::ENodeType::INPUT) {
				throw std::runtime_error("Input node cannot be output of connection");
			}
		}
	}
}

// Validate network topology (no cycles in feedforward networks)
void NetworkValidator::validateTopology(const core::Genome& genome) {
	std::map<int, std::vector<int>> adjacencyList;
	std::set<int> visited;
	std::set<int> recursionStack;
	
	// Build adjacency list from enabled genes
	for (const auto& gene : genome.getGenes()) {
		if (gene.enabled) {
			adjacencyList[gene.inputNode].push_back(gene.outputNode);
		}
	}
	
	// Check for cycles using DFS
	std::function<bool(int)> hasCycle = [&](int node) {
		visited.insert(node);
		recursionStack.insert(node);
		
		for (int neighbor : adjacencyList[node]) {
			if (visited.find(neighbor) == visited.end()) {
				if (hasCycle(neighbor)) return true;
			} else if (recursionStack.find(neighbor) != recursionStack.end()) {
				return true;
			}
		}
		
		recursionStack.erase(node);
		return false;
	};
	
	for (const auto& [nodeId, _] : genome.getNodes()) {
		if (visited.find(nodeId) == visited.end()) {
			if (hasCycle(nodeId)) {
				throw std::runtime_error("Cycle detected in feedforward network");
			}
		}
	}
}

// Validate all output nodes are reachable from inputs
void NetworkValidator::validateConnectivity(const core::Genome& genome) {
	std::set<int> reachable;
	std::queue<int> queue;
	
	// Start with inputs and bias
	for (const auto& [id, type] : genome.getNodes()) {
		if (type == core::ENodeType::INPUT || id == -1) {
			queue.push(id);
			reachable.insert(id);
		}
	}
	
	// Breadth-first search through enabled connections
	while (!queue.empty()) {
		int current = queue.front();
		queue.pop();
		
		for (const auto& gene : genome.getGenes()) {
			if (gene.enabled && gene.inputNode == current) {
				if (reachable.insert(gene.outputNode).second) {
					queue.push(gene.outputNode);
				}
			}
		}
	}
	
	// Verify all output nodes are reachable
	for (const auto& [id, type] : genome.getNodes()) {
		if (type == core::ENodeType::OUTPUT && reachable.find(id) == reachable.end()) {
			throw std::runtime_error("Output node " + std::to_string(id) + " is not reachable");
		}
	}
}

}
}
