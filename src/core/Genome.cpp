#include "Genome.hpp"
#include "Gene.hpp"
#include "InnovationTracker.hpp"
#include <random>
#include <algorithm>
#include <stdexcept>

namespace neat {
namespace core {

Genome::Genome(const Config& config)
    : config(config)
    , network(std::make_unique<network::Network>(config.networkConfig))
    , rng(std::random_device{}()) {}

void Genome::initMinimalTopology(int32_t inputs, int32_t outputs) {
    // Add bias node
    addNode(-1, ENodeType::BIAS);
    
    // Add input and output nodes
    for (int32_t i = 0; i < inputs; ++i) {
        addNode(getNextNodeId(), ENodeType::INPUT);
    }
    
    for (int32_t i = 0; i < outputs; ++i) {
        addNode(getNextNodeId(), ENodeType::OUTPUT);
    }
    
    // Connect inputs to outputs
    std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
    
    for (const auto& [fromId, fromType] : nodes) {
        if (fromType != ENodeType::INPUT && fromType != ENodeType::BIAS) continue;
        
        for (const auto& [toId, toType] : nodes) {
            if (toType != ENodeType::OUTPUT) continue;
            addConnection(fromId, toId, weightDist(rng));
        }
    }
}

void Genome::addNode(NodeId id, ENodeType type) {
    nodes[id] = type;
    maxNodeId = std::max(maxNodeId, id);
    network->addNode(id, type);
}

void Genome::addGene(const Gene& gene) {
    genes.push_back(gene);
    network->addConnection(gene.inputNode, gene.outputNode, gene.weight, gene.enabled);
}

void Genome::addConnection(NodeId from, NodeId to, double weight) {
    Gene gene(from, to, weight, true, InnovationTracker::getNextInnovation());
    addGene(gene);
}

bool Genome::addConnectionMutation() {
    std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
    
    // Find all possible connections
    std::vector<std::pair<NodeId, NodeId>> possibleConnections;
    
    for (const auto& [fromId, fromType] : nodes) {
        if (fromType == ENodeType::OUTPUT) continue;
        
        for (const auto& [toId, toType] : nodes) {
            if (toType == ENodeType::INPUT || fromId == toId) continue;
            
            if (isValidConnection(fromId, toId)) {
                possibleConnections.emplace_back(fromId, toId);
            }
        }
    }
    
    if (possibleConnections.empty()) return false;
    
    // Select random connection
    std::uniform_int_distribution<size_t> connDist(0, possibleConnections.size() - 1);
    auto [fromId, toId] = possibleConnections[connDist(rng)];
    
    addConnection(fromId, toId, weightDist(rng));
    return true;
}

bool Genome::addNodeMutation() {
    if (genes.empty()) return false;
    
    // Select random enabled gene
    std::vector<std::reference_wrapper<Gene>> enabledGenes;
    for (auto& gene : genes) {
        if (gene.enabled) {
            enabledGenes.push_back(gene);
        }
    }
    
    if (enabledGenes.empty()) return false;
    
    std::uniform_int_distribution<size_t> geneDist(0, enabledGenes.size() - 1);
    Gene& selectedGene = enabledGenes[geneDist(rng)];
    
    // Disable selected gene
    selectedGene.enabled = false;
    
    // Add new node
    NodeId newNodeId = getNextNodeId();
    addNode(newNodeId, ENodeType::HIDDEN);
    
    // Add new connections
    addConnection(selectedGene.inputNode, newNodeId, 1.0);
    addConnection(newNodeId, selectedGene.outputNode, selectedGene.weight);
    
    return true;
}

void Genome::mutateWeights() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::normal_distribution<double> perturbDist(0.0, config.weightPerturbationRange);
    std::uniform_real_distribution<double> newWeightDist(-config.newWeightRange, config.newWeightRange);
    
    for (auto& gene : genes) {
        if (dist(rng) < config.weightMutationRate) {
            if (dist(rng) < config.weightPerturbationRate) {
                gene.weight += perturbDist(rng);
            } else {
                gene.weight = newWeightDist(rng);
            }
        }
    }
}

Genome Genome::crossover(const Genome& parent1, const Genome& parent2, const Config& config) {
    Genome child(config);
    
    // Inherit nodes
    child.nodes = parent1.nodes;
    
    // Match genes by innovation number
    for (const auto& gene1 : parent1.genes) {
        auto it = std::find_if(parent2.genes.begin(), parent2.genes.end(),
            [&gene1](const Gene& gene2) { return gene1.innovation == gene2.innovation; });
            
        if (it != parent2.genes.end()) {
            // Matching gene - randomly inherit from either parent
            child.addGene(std::uniform_int_distribution<int>(0,1)(child.rng) ? gene1 : *it);
        } else {
            // Disjoint/excess gene - inherit from more fit parent
            if (parent1.fitness >= parent2.fitness) {
                child.addGene(gene1);
            }
        }
    }
    
    return child;
}

double Genome::compatibilityDistance(const Genome& genome1, const Genome& genome2) {
    double disjoint = 0.0;
    double weightDiff = 0.0;
    int matching = 0;
    
    size_t i = 0, j = 0;
    while (i < genome1.genes.size() && j < genome2.genes.size()) {
        if (genome1.genes[i].innovation == genome2.genes[j].innovation) {
            weightDiff += std::abs(genome1.genes[i].weight - genome2.genes[j].weight);
            matching++;
            i++;
            j++;
        } else if (genome1.genes[i].innovation < genome2.genes[j].innovation) {
            disjoint++;
            i++;
        } else {
            disjoint++;
            j++;
        }
    }
    
    disjoint += (genome1.genes.size() - i) + (genome2.genes.size() - j);
    weightDiff = matching > 0 ? weightDiff / matching : 0;
    
    const double c1 = 1.0, c2 = 1.0;
    return (c1 * disjoint / std::max(genome1.genes.size(), genome2.genes.size())) + 
           (c2 * weightDiff);
}

bool Genome::isValidConnection(NodeId from, NodeId to) const {
    // Check if connection already exists
    auto it = std::find_if(genes.begin(), genes.end(),
        [from, to](const Gene& gene) {
            return gene.inputNode == from && gene.outputNode == to;
        });
    
    return it == genes.end();
}

void Genome::rebuildNetwork() {
    network = std::make_unique<network::Network>(config.networkConfig);
    
    // Rebuild nodes
    for (const auto& [id, type] : nodes) {
        network->addNode(id, type);
    }
    
    // Rebuild connections
    for (const auto& gene : genes) {
        network->addConnection(gene.inputNode, gene.outputNode, gene.weight, gene.enabled);
    }
}

std::vector<double> Genome::activate(const std::vector<double>& inputs) const {
    return network->activate(inputs);
}

bool Genome::validate() const {
    return network->validate();
}

}
}
