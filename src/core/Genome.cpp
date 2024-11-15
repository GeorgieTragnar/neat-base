#include "Genome.hpp"
#include "Gene.hpp"
#include "InnovationTracker.hpp"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <iostream>

namespace neat {
namespace core {


// Validation helper methods
bool Genome::validateNodeStructure() const {
    // Count input and output nodes
    int inputCount = 0;
    int outputCount = 0;
    bool hasBias = false;
    
    for (const auto& [id, type] : nodes) {
        if (id == -1 && type == ENodeType::BIAS) {
            hasBias = true;
        } else if (type == ENodeType::INPUT) {
            inputCount++;
        } else if (type == ENodeType::OUTPUT) {
            outputCount++;
        }
    }
    
    // Network must have:
    // 1. At least one input node
    // 2. At least one output node
    // 3. A bias node
    return inputCount > 0 && outputCount > 0 && hasBias;
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
    std::cout << "Initializing minimal topology with " << inputs 
              << " inputs and " << outputs << " outputs" << std::endl;
              
    // Clear any existing structure
    nodes.clear();
    genes.clear();
    maxNodeIdx = -1;
    
    addNode(-1, ENodeType::BIAS, false);
    std::cout << "Added bias node" << std::endl;

    // Add input nodes
    std::vector<NodeId> inputNodeIds;
    for (int32_t i = 0; i < inputs; ++i) {
        NodeId id = getNextNodeId();
        addNode(id, ENodeType::INPUT, false);
        inputNodeIds.push_back(id);
    }
    std::cout << "Added " << inputs << " input nodes" << std::endl;
    
    // Add output nodes
    std::vector<NodeId> outputNodeIds;
    for (int32_t i = 0; i < outputs; ++i) {
        NodeId id = getNextNodeId();
        addNode(id, ENodeType::OUTPUT, false);
        outputNodeIds.push_back(id);
    }
    std::cout << "Added " << outputs << " output nodes" << std::endl;
    
    // Add connections from all inputs (including bias) to all outputs
    std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
    
    // Connect regular inputs to outputs
    for (NodeId inputId : inputNodeIds) {
        for (NodeId outputId : outputNodeIds) {
            addConnection(inputId, outputId, weightDist(rng), false);
        }
    }

    // Connect bias to all outputs
    for (NodeId outputId : outputNodeIds) {
        addConnection(-1, outputId, weightDist(rng), false);
    }

    std::cout << "Added connections between nodes" << std::endl;
    
    // Print final network structure
    std::cout << "Final network structure:" << std::endl;
    std::cout << "Nodes:" << std::endl;
    for (const auto& [id, type] : nodes) {
        std::cout << "  Node " << id << ": " 
                  << (type == ENodeType::INPUT ? "INPUT" :
                      type == ENodeType::OUTPUT ? "OUTPUT" :
                      type == ENodeType::BIAS ? "BIAS" : "HIDDEN")
                  << std::endl;
    }

    // Now validate the complete network
    rebuildNetwork();
    validate();
}

// Add these validation helper methods to the Genome class private section
void Genome::ensureNodeExists(NodeId id, ENodeType type) {
    if (nodes.find(id) == nodes.end()) {
        nodes[id] = type;
        if (id > maxNodeIdx && id != -1) {
            maxNodeIdx = id;
        }
        // Update network
        network->addNode(id, type);
    }
}

std::vector<std::pair<Genome::NodeId, Genome::NodeId>> Genome::findPossibleConnections() const {
    std::vector<std::pair<NodeId, NodeId>> result;
    
    // Create set of existing connections for quick lookup
    std::set<std::pair<NodeId, NodeId>> existingConns;
    for (const auto& gene : genes) {
        if (gene.enabled) {
            existingConns.insert({gene.inputNode, gene.outputNode});
        }
    }
    
    // Check all possible node pairs
    for (const auto& [fromId, fromType] : nodes) {
        if (fromType == ENodeType::OUTPUT) continue; // Outputs can't be source nodes
        
        for (const auto& [toId, toType] : nodes) {
            if (toType == ENodeType::INPUT) continue; // Inputs can't be target nodes
            if (fromId == toId) continue; // No self-connections
            
            // Skip if connection already exists
            if (existingConns.find({fromId, toId}) != existingConns.end()) {
                continue;
            }
            
            // Skip if would create cycle (for feed-forward networks)
            if (!config.networkConfig.allowRecurrent && wouldCreateCycle(fromId, toId)) {
                continue;
            }
            
            result.emplace_back(fromId, toId);
        }
    }
    
    return result;
}

// Add this implementation to Genome.cpp:
bool Genome::wouldCreateCycle(NodeId fromId, NodeId toId) const {
    if (fromId == toId) {
        return true;  // Self-connection is a cycle
    }
    
    // Create adjacency list for enabled connections
    std::map<NodeId, std::vector<NodeId>> adjacencyList;
    for (const auto& gene : genes) {
        if (gene.enabled) {
            adjacencyList[gene.inputNode].push_back(gene.outputNode);
        }
    }
    
    // Temporarily add the new connection
    adjacencyList[fromId].push_back(toId);
    
    // Use DFS to detect cycles
    std::set<NodeId> visited;
    std::set<NodeId> recursionStack;
    
    std::function<bool(NodeId)> hasCycleUtil = [&](NodeId node) -> bool {
        if (recursionStack.find(node) != recursionStack.end()) {
            return true;  // Node is in recursion stack - cycle found
        }
        
        if (visited.find(node) != visited.end()) {
            return false;  // Already visited this path - no cycle here
        }
        
        visited.insert(node);
        recursionStack.insert(node);
        
        // Check all neighbors
        for (NodeId neighbor : adjacencyList[node]) {
            if (hasCycleUtil(neighbor)) {
                return true;
            }
        }
        
        recursionStack.erase(node);
        return false;
    };
    
    // Start DFS from the new source node
    return hasCycleUtil(fromId);
}

void Genome::addNode(NodeId id, ENodeType type, bool validateAfter) {
    // Don't allow duplicate nodes
    if (nodes.find(id) != nodes.end()) {
        return;
    }
    
    nodes[id] = type;
    if (id > maxNodeIdx && id != -1) {
        maxNodeIdx = id;
    }
    
    // Keep network in sync
    network->addNode(id, type);
    
    if (validateAfter) {
        validate();
    }
}

void Genome::addGene(const Gene& gene, bool validateAfter) {
    // Check if gene already exists
    auto it = std::find_if(genes.begin(), genes.end(),
        [&gene](const Gene& existing) {
            return existing.inputNode == gene.inputNode && 
                   existing.outputNode == gene.outputNode;
        });
        
    if (it != genes.end()) {
        // Update existing gene
        it->weight = gene.weight;
        it->enabled = gene.enabled;
        it->innovation = gene.innovation;
        
        // Update connection in network if enabled status or weight changed
        if (gene.enabled) {
            network->addConnection(gene.inputNode, gene.outputNode, gene.weight);
        }
    } else {
        // Add new gene
        genes.push_back(gene);
        
        // Add to network if enabled
        if (gene.enabled) {
            network->addConnection(gene.inputNode, gene.outputNode, gene.weight);
        }
    }
    
    if (validateAfter) {
        validate();
    }
}

void Genome::addConnection(NodeId from, NodeId to, double weight, bool validateAfter) {
    // Ensure both nodes exist before adding connection
    if (from != -1) {  // Skip check for bias node
        auto fromType = nodes.find(from) != nodes.end() ? nodes[from] : ENodeType::HIDDEN;
        ensureNodeExists(from, fromType);
    }
    
    auto toType = nodes.find(to) != nodes.end() ? nodes[to] : ENodeType::HIDDEN;
    ensureNodeExists(to, toType);
    
    // Now add the connection
    Gene gene(from, to, weight, true, InnovationTracker::getNextInnovation());
    genes.push_back(gene);
    
    // Update network
    if (gene.enabled) {
        network->addConnection(gene.inputNode, gene.outputNode, gene.weight);
    }
    
    if (validateAfter) {
        validate();
    }
}

bool Genome::addConnectionMutation() {
    std::cout << "Starting connection mutation with " 
              << std::count_if(nodes.begin(), nodes.end(),
                  [](const auto& pair) { return pair.second == ENodeType::INPUT; })
              << " inputs and " 
              << std::count_if(nodes.begin(), nodes.end(),
                  [](const auto& pair) { return pair.second == ENodeType::OUTPUT; })
              << " outputs" << std::endl;
              
    // Make a backup before mutation
    auto backup = *this;
    
    try {
        auto possibleConnections = findPossibleConnections();
        if (possibleConnections.empty()) {
            return false;
        }
        
        // Select random possible connection
        std::uniform_int_distribution<size_t> dist(0, possibleConnections.size() - 1);
        auto [fromId, toId] = possibleConnections[dist(rng)];
        
        // Generate random weight
        std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
        double weight = weightDist(rng);
        
        // Add new connection with validation
        addConnection(fromId, toId, weight);
        
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
    
    // Debug output for node structure
    std::cout << "Validating genome..." << std::endl;
    std::cout << "Node structure:" << std::endl;
    for (const auto& [id, type] : nodes) {
        std::cout << "  Node " << id << ": " 
                  << (type == ENodeType::INPUT ? "INPUT" :
                      type == ENodeType::OUTPUT ? "OUTPUT" :
                      type == ENodeType::BIAS ? "BIAS" : "HIDDEN")
                  << std::endl;
    }

    // Check node structure
    if (!validateNodeStructure()) {
        int inputCount = 0;
        int outputCount = 0;
        int biasCount = 0;
        for (const auto& [id, type] : nodes) {
            if (type == ENodeType::INPUT) inputCount++;
            if (type == ENodeType::OUTPUT) outputCount++;
            if (type == ENodeType::BIAS) biasCount++;
        }
        errors.push_back("Invalid node structure: Found " + 
                        std::to_string(inputCount) + " inputs, " +
                        std::to_string(outputCount) + " outputs, and " +
                        std::to_string(biasCount) + " bias nodes");
    }
    
    // Debug output for genes
    std::cout << "Gene structure:" << std::endl;
    for (const auto& gene : genes) {
        std::cout << "  Gene: " << gene.inputNode << " -> " << gene.outputNode 
                  << " (weight: " << gene.weight 
                  << ", enabled: " << (gene.enabled ? "yes" : "no") << ")" << std::endl;
    }
    
    // Check gene structure
    if (!validateGeneStructure()) {
        errors.push_back("Invalid gene structure");
    }
    
    // Check for cycles in feedforward networks
    if (!config.networkConfig.allowRecurrent && hasCycle()) {
        errors.push_back("Network contains cycles in feedforward mode");
    }
    
    // Check output node reachability
    auto reachable = getReachableNodes();
    std::cout << "Reachable nodes: ";
    for (const auto& id : reachable) {
        std::cout << id << " ";
    }
    std::cout << std::endl;
    
    for (const auto& [id, type] : nodes) {
        if (type == ENodeType::OUTPUT && reachable.find(id) == reachable.end()) {
            errors.push_back("Output node " + std::to_string(id) + " is not reachable");
        }
    }
    
    // Network level validation
    if (network && !network->validate()) {
        errors.push_back("Neural network validation failed");
    }
    
    if (!errors.empty()) {
        std::string errorMsg = "Network validation failed:\n";
        for (const auto& error : errors) {
            errorMsg += "- " + error + "\n";
        }
        std::cout << "Validation failed with errors:\n" << errorMsg;
        throw std::runtime_error(errorMsg);
    }
    
    std::cout << "Validation successful" << std::endl;
    return true;
}

}
}
