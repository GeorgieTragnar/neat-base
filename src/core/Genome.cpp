#include "Genome.hpp"
#include "Gene.hpp"
#include "InnovationTracker.hpp"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <queue>

namespace neat {
namespace core {


// Validation helper methods
bool Genome::validateNodeStructure() const {
    // Check node ID bounds and required types
    bool hasInput = false, hasOutput = false;
    for (const auto& [id, type] : nodes) {
        if (id > maxNodeIdx && id != -1) {
            return false;
        }
        if (type == ENodeType::INPUT) hasInput = true;
        if (type == ENodeType::OUTPUT) hasOutput = true;
    }
    return hasInput && hasOutput;
}

bool Genome::validateGeneStructure() const {
    for (const auto& gene : genes) {
        if (gene.enabled) {
            // Validate node existence
            if (gene.inputNode != -1 && nodes.find(gene.inputNode) == nodes.end()) {
                return false;
            }
            if (nodes.find(gene.outputNode) == nodes.end()) {
                return false;
            }

            // Validate node types
            if (gene.inputNode != -1) {
                auto inputType = nodes.at(gene.inputNode);
                if (inputType == ENodeType::OUTPUT) return false;
            }
            auto outputType = nodes.at(gene.outputNode);
            if (outputType == ENodeType::INPUT) return false;
        }
    }
    return true;
}

bool Genome::hasCycle() const {
    std::map<int32_t, std::vector<int32_t>> adjacencyList;
    std::set<int32_t> visited;
    std::set<int32_t> recursionStack;
    
    // Build adjacency list from enabled genes
    for (const auto& gene : genes) {
        if (gene.enabled) {
            adjacencyList[gene.inputNode].push_back(gene.outputNode);
        }
    }
    
    std::function<bool(int32_t)> hasCycleUtil = [&](int32_t node) {
        visited.insert(node);
        recursionStack.insert(node);
        
        for (int32_t neighbor : adjacencyList[node]) {
            if (visited.find(neighbor) == visited.end()) {
                if (hasCycleUtil(neighbor)) return true;
            } else if (recursionStack.find(neighbor) != recursionStack.end()) {
                return true;
            }
        }
        
        recursionStack.erase(node);
        return false;
    };
    
    for (const auto& [nodeId, _] : nodes) {
        if (visited.find(nodeId) == visited.end()) {
            if (hasCycleUtil(nodeId)) return true;
        }
    }
    
    return false;
}

std::set<int32_t> Genome::getReachableNodes() const {
    std::set<int32_t> reachable;
    std::queue<int32_t> queue;
    
    // Start with inputs and bias
    for (const auto& [id, type] : nodes) {
        if (type == ENodeType::INPUT || id == -1) {
            queue.push(id);
            reachable.insert(id);
        }
    }
    
    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();
        
        for (const auto& gene : genes) {
            if (gene.enabled && gene.inputNode == current) {
                if (reachable.insert(gene.outputNode).second) {
                    queue.push(gene.outputNode);
                }
            }
        }
    }
    return reachable;
}

void Genome::sanitizeNetwork() {
    // First disable invalid genes
    for (auto& gene : genes) {
        if (gene.enabled) {
            if ((gene.inputNode != -1 && nodes.find(gene.inputNode) == nodes.end()) ||
                nodes.find(gene.outputNode) == nodes.end()) {
                gene.enabled = false;
            }
        }
    }

    // Remove cycles
    while (hasCycle()) {
        for (auto it = genes.rbegin(); it != genes.rend(); ++it) {
            if (it->enabled) {
                it->enabled = false;
                if (!hasCycle()) break;
            }
        }
    }

    // Ensure connectivity
    auto reachable = getReachableNodes();
    for (auto& gene : genes) {
        if (gene.enabled &&
            (reachable.find(gene.inputNode) == reachable.end() ||
                reachable.find(gene.outputNode) == reachable.end())) {
            gene.enabled = false;
        }
    }
}

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

void Genome::addNode(NodeId id, ENodeType type, bool validateAfter) {
    nodes[id] = type;
    if (id > maxNodeIdx && id != -1) {
        maxNodeIdx = id;
    }
    network->addNode(id, type);
    if (validateAfter)
        validate();
}

void Genome::addGene(const Gene& gene, bool validateAfter) {
    genes.push_back(gene);
    if (gene.enabled)
        network->addConnection(gene.inputNode, gene.outputNode, gene.weight, gene.enabled);
    if (validateAfter)
        validate();
}

void Genome::addConnection(NodeId from, NodeId to, double weight) {
    Gene gene(from, to, weight, true, InnovationTracker::getNextInnovation());
    addGene(gene);
}

bool Genome::addConnectionMutation() {
    // Make a backup before mutation
    auto backup = *this;
    
    try {
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
        
        // Validate the modified network
        validate();
        return true;
        
    } catch (const std::exception& e) {
        // Restore backup on failure
        *this = backup;
        return false;
    }
}

// Clone method for safe copying
Genome Genome::clone() const {
    Genome copy(*this);
    copy.validate();
    return copy;
}

bool Genome::addNodeMutation() {
    // Make a backup before mutation
    auto backup = *this;
    
    try {
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
        
        // Validate the modified network
        validate();
        return true;
        
    } catch (const std::exception& e) {
        // Restore backup on failure
        *this = backup;
        return false;
    }
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

// Modify rebuildNetwork to ensure proper initialization
void Genome::rebuildNetwork() {
    network = std::make_unique<network::Network>(config.networkConfig);
    
    // First add all nodes
    for (const auto& [id, type] : nodes) {
        network->addNode(id, type);
    }
    
    // Then add all connections
    for (const auto& gene : genes) {
        if (gene.enabled) {
            network->addConnection(gene.inputNode, gene.outputNode, gene.weight);
        }
    }
}

std::vector<double> Genome::activate(const std::vector<double>& inputs) const {
    validate();
    return network->activate(inputs);
}

bool Genome::validate() const {
    std::vector<std::string> errors;
    
    if (!validateNodeStructure()) {
        errors.push_back("Invalid node structure");
    }
    if (!validateGeneStructure()) {
        errors.push_back("Invalid gene structure");
    }
    if (hasCycle()) {
        errors.push_back("Network contains cycles");
    }
    
    auto reachable = getReachableNodes();
    for (const auto& [id, type] : nodes) {
        if (type == ENodeType::OUTPUT && reachable.find(id) == reachable.end()) {
            errors.push_back("Output node " + std::to_string(id) + " is not reachable");
        }
    }
    
    // Add network validation
    if (!network->validate()) {
        errors.push_back("Neural network validation failed");
    }
    
    if (!errors.empty()) {
        std::string errorMsg = "Network validation failed:\n";
        for (const auto& error : errors) {
            errorMsg += "- " + error + "\n";
        }
        throw std::runtime_error(errorMsg);
    }

    return network->validate();
}

}
}
