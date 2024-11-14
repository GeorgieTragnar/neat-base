
#include "MutationOperator.hpp"
#include "core/Gene.hpp"

namespace neat {
namespace evolution {

MutationOperator::MutationOperator(const Config& config)
        : config(config)
        , rng(std::random_device{}())
        , weightDist(0.0, 1.0)
        , perturbDist(0.0, config.weightPerturbationRange)
        , newWeightDist(0.0, config.newWeightRange)
    {}

void MutationOperator::mutate(core::Genome& genome) {
	// Weight mutations
	if (weightDist(rng) < config.weightMutationRate) {
		for (auto& gene : genome.getGenes()) {
			if (weightDist(rng) < config.weightPerturbationRate) {
				// Perturb existing weight
				gene.weight += perturbDist(rng);
			} else {
				// Assign new random weight
				gene.weight = newWeightDist(rng);
			}
		}
	}

	// Structural mutations
	if (weightDist(rng) < config.newNodeRate) {
		addNodeMutation(genome);
	}

	if (weightDist(rng) < config.newConnectionRate) {
		addConnectionMutation(genome);
	}
}

bool MutationOperator::addNodeMutation(core::Genome& genome) {
	auto& genes = genome.getGenes();
	if (genes.empty()) return false;

	// Get enabled genes
	std::vector<std::reference_wrapper<core::Gene>> enabledGenes;
	for (auto& gene : genes) {
		if (gene.enabled) {
			enabledGenes.push_back(gene);
		}
	}
	
	if (enabledGenes.empty()) return false;

	// Select random enabled gene
	std::uniform_int_distribution<size_t> geneDist(0, enabledGenes.size() - 1);
	core::Gene& selectedGene = enabledGenes[geneDist(rng)];
	
	// Create new node
	int32_t newNodeId = genome.getNextNodeId();
	double oldWeight = selectedGene.weight;
	
	// Disable old connection
	selectedGene.enabled = false;
	
	// Add new node and connections
	genome.addNode(newNodeId, core::ENodeType::HIDDEN);
	genome.addConnection(selectedGene.inputNode, newNodeId, 1.0);
	genome.addConnection(newNodeId, selectedGene.outputNode, oldWeight);
	
	return true;
}

bool MutationOperator::addConnectionMutation(core::Genome& genome) {
	auto possibleConnections = findPossibleConnections(genome);
	if (possibleConnections.empty()) return false;

	// Select random possible connection
	std::uniform_int_distribution<size_t> connDist(0, possibleConnections.size() - 1);
	auto [fromNode, toNode] = possibleConnections[connDist(rng)];
	
	// Create new connection with random weight
	return genome.addConnection(fromNode, toNode, newWeightDist(rng));
}

std::vector<std::pair<int32_t, int32_t>> MutationOperator::findPossibleConnections(const core::Genome& genome) {
	std::vector<std::pair<int32_t, int32_t>> connections;
	const auto& nodes = genome.getNodes();

	for (const auto& [fromId, fromType] : nodes) {
		if (fromType == core::ENodeType::OUTPUT) continue;

		for (const auto& [toId, toType] : nodes) {
			if (toType == core::ENodeType::INPUT || fromId == toId) continue;
			
			// Check if connection would be valid
			if (genome.isValidConnection(fromId, toId)) {
				connections.emplace_back(fromId, toId);
			}
		}
	}

	// Also consider bias node connections
	for (const auto& [nodeId, nodeType] : nodes) {
		if (nodeType != core::ENodeType::INPUT && genome.isValidConnection(-1, nodeId)) {
			connections.emplace_back(-1, nodeId);
		}
	}

	return connections;
}

}
}
