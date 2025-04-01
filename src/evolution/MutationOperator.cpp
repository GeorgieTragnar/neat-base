
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

	// Activation function mutations
	if (weightDist(rng) < config.activationMutationRate) {
		mutateActivations(genome);
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
    
    // Create new activation genes for the split connections
    auto inputActivation = selectedGene.activation; // Keep original activation for input
    auto outputActivation = core::ActivationGene::createRandom(config.activationConfig);
    
    // Add new node and connections
    genome.addNode(newNodeId, core::ENodeType::HIDDEN, true);
    
    // Add two new connections with appropriate activations
    genome.addConnection(
        selectedGene.inputNode, 
        newNodeId, 
        1.0,  // Weight to new node
        true, 
        inputActivation.getType()
    );
    
    genome.addConnection(
        newNodeId,
        selectedGene.outputNode,
        oldWeight,  // Original weight to output
        true,
        outputActivation.getType()
    );

	return true;
}

bool MutationOperator::addConnectionMutation(core::Genome& genome) {
	auto possibleConnections = findPossibleConnections(genome);
	if (possibleConnections.empty()) return false;

	// Select random possible connection
	std::uniform_int_distribution<size_t> connDist(0, possibleConnections.size() - 1);
    auto [fromId, toId] = possibleConnections[connDist(rng)];
	
    // Generate random weight
    std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
    double weight = weightDist(rng);

	// Create new activation function based on node types
	const auto& nodes = genome.getNodes();
	auto toType = nodes.at(toId);
	
	core::ActivationGene activation;
	if (toType == core::ENodeType::OUTPUT) {
		// For output nodes, use only basic tier activations
		std::uniform_int_distribution<int> basicDist(0, 2); // SIGMOID, TANH, RELU
		activation = core::ActivationGene(static_cast<core::EActivationType>(basicDist(rng)));
	} else {
		// For hidden nodes, allow any activation function
		activation = core::ActivationGene::createRandom(config.activationConfig);
	}

	// Create new connection with random weight
	genome.addConnection(fromId, toId, weight, true, activation.getType());
	
	return true;
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


// MutationOperator.cpp - Implementation of activation mutations
void MutationOperator::mutateActivations(core::Genome& genome) {
    for (auto& gene : genome.getGenes()) {
        if (!gene.enabled) continue;
        
        if (weightDist(rng) < config.activationMutationRate) {
            mutateConnectionActivation(gene, genome);
        }
    }
}

void MutationOperator::mutateConnectionActivation(core::Gene& gene, core::Genome& genome) {
    const auto& nodes = genome.getNodes();
    auto targetType = nodes.at(gene.outputNode);
    auto currentType = gene.activation.getType();
    
    // Special handling for output nodes - restrict to basic tier
    if (targetType == core::ENodeType::OUTPUT) {
        auto newType = selectNewActivation(currentType, targetType);
        gene.activation = core::ActivationGene(newType);
        return;
    }
    
    int currentTier = getCurrentTierIndex(currentType);
    
    // Determine if we should change tiers
    if (weightDist(rng) < config.interTierMutationRate) {
        if (shouldMutateToHigherTier(gene.activation)) {
            currentTier = std::min(currentTier + 1, static_cast<int>(tiers.size()) - 1);
        } else if (shouldMutateToLowerTier(gene.activation)) {
            currentTier = std::max(currentTier - 1, 0);
        }
    }
    
    // Select new activation from the determined tier
    const auto& tierFunctions = tiers[currentTier].functions;
    std::uniform_int_distribution<size_t> dist(0, tierFunctions.size() - 1);
    auto newType = tierFunctions[dist(rng)];
    
    // Avoid selecting the same function
    if (newType == currentType && tierFunctions.size() > 1) {
        do {
            newType = tierFunctions[dist(rng)];
        } while (newType == currentType);
    }
    
    gene.activation = core::ActivationGene(newType);
}

core::EActivationType MutationOperator::selectNewActivation(
    core::EActivationType current, 
    core::ENodeType targetType) {
    
    // For output nodes, restrict to basic tier
    if (targetType == core::ENodeType::OUTPUT) {
        const auto& basicFunctions = tiers[0].functions;
        std::uniform_int_distribution<size_t> dist(0, basicFunctions.size() - 1);
        
        auto newType = basicFunctions[dist(rng)];
        if (newType == current && basicFunctions.size() > 1) {
            do {
                newType = basicFunctions[dist(rng)];
            } while (newType == current);
        }
        return newType;
    }
    
    // For hidden nodes, use tier-based selection
    double prob = weightDist(rng);
    double cumProb = 0.0;
    
    for (const auto& tier : tiers) {
        cumProb += tier.mutationProbability;
        if (prob <= cumProb) {
            std::uniform_int_distribution<size_t> dist(0, tier.functions.size() - 1);
            return tier.functions[dist(rng)];
        }
    }
    
    // Fallback to basic tier
    return tiers[0].functions[0];
}

bool MutationOperator::shouldMutateToHigherTier(const core::ActivationGene& current) {
    int currentTier = getCurrentTierIndex(current.getType());
    
    // Only allow moving up if we're not in the highest tier
    if (currentTier >= static_cast<int>(tiers.size()) - 1) {
        return false;
    }
    
    // Higher tiers require better fitness or specific conditions
    double advancementProb = config.activationConfig.interTierMutationRate;
    
    // Adjust probability based on current tier
    if (currentTier == 0) {
        advancementProb *= 0.5;  // Harder to move from basic tier
    }
    
    return weightDist(rng) < advancementProb;
}

bool MutationOperator::shouldMutateToLowerTier(const core::ActivationGene& current) {
    int currentTier = getCurrentTierIndex(current.getType());
    
    // Only allow moving down if we're not in the lowest tier
    if (currentTier <= 0) {
        return false;
    }
    
    // Higher chance to move down from experimental tier
    double regressionProb = config.activationConfig.interTierMutationRate;
    if (currentTier == 2) {
        regressionProb *= 1.5;  // More likely to move down from experimental
    }
    
    return weightDist(rng) < regressionProb;
}

int MutationOperator::getCurrentTierIndex(core::EActivationType type) {
    for (size_t i = 0; i < tiers.size(); ++i) {
        const auto& functions = tiers[i].functions;
        if (std::find(functions.begin(), functions.end(), type) != functions.end()) {
            return i;
        }
    }
    return 0;  // Default to basic tier
}

}
}
