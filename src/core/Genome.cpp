#include "Genome.hpp"
#include "Gene.hpp"
#include "InnovationTracker.hpp"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <iostream>

#include "logger/Logger.hpp"
static auto logger = LOGGER("core::Genome");

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
    LOG_DEBUG("Initializing minimal topology with {} inputs and {} outputs", inputs, outputs);
              
    // Clear any existing structure
    nodes.clear();
    genes.clear();
    maxNodeIdx = -1;
    
    addNode(-1, ENodeType::BIAS, false);
    LOG_TRACE("Added bias node");
    
    // Add input nodes
    std::vector<NodeId> inputNodeIds;
    for (int32_t i = 0; i < inputs; ++i) {
        NodeId id = getNextNodeId();
        addNode(id, ENodeType::INPUT, false);
        inputNodeIds.push_back(id);
    }
    LOG_TRACE("Added {} input nodes", inputs);
    
    // Add output nodes
    std::vector<NodeId> outputNodeIds;
    for (int32_t i = 0; i < outputs; ++i) {
        NodeId id = getNextNodeId();
        addNode(id, ENodeType::OUTPUT, false);
        outputNodeIds.push_back(id);
    }
    LOG_TRACE("Added {} output nodes", outputs);
    
    // Add hidden nodes
    std::vector<NodeId> hiddenNodeIds;
    for (int32_t i = 0; i < config.maxHiddenNodes; ++i) {
        NodeId id = getNextNodeId();
        addNode(id, ENodeType::HIDDEN, false);
        hiddenNodeIds.push_back(id);
    }
    LOG_TRACE("Added {} hidden nodes", config.maxHiddenNodes);

    // Initialize with minimal connections - one random input per output
    std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
    
    // Initialize connections through hidden layer
    for (NodeId hiddenId : hiddenNodeIds) {
        // Connect inputs to hidden
        for (NodeId inputId : inputNodeIds) {
            auto activation = core::ActivationGene::createRandom(config.activationConfig);
            addConnection(inputId, hiddenId, weightDist(rng), false, activation.getType());
        }
        auto biasActivation = core::ActivationGene::createRandom(config.activationConfig);
        addConnection(-1, hiddenId, weightDist(rng), false, biasActivation.getType()); // bias to hidden
        
        // Connect hidden to outputs
        for (NodeId outputId : outputNodeIds) {
            auto activation = core::ActivationGene::createRandom(config.activationConfig);
            addConnection(hiddenId, outputId, weightDist(rng), false, activation.getType());
        }
    }
    
    for (NodeId outputId : outputNodeIds) {
        // Connect regular inputs
        for (NodeId inputId : inputNodeIds) {
            auto activation = core::ActivationGene::createRandom(config.activationConfig);
            addConnection(inputId, outputId, weightDist(rng), false, activation.getType());
            LOG_TRACE("Added initial connection from {} to {}", inputId, outputId);
        }
        // Connect bias
        auto biasActivation = core::ActivationGene::createRandom(config.activationConfig);
        addConnection(-1, outputId, weightDist(rng), false, biasActivation.getType());
        LOG_TRACE("Added initial connection from {} to {}", -1, outputId);
    }

    LOG_TRACE("Added minimal initial connections");
    
    // Print final network structure
    LOG_INFO("Final network structure:");
    LOG_INFO("Nodes:");
    for (const auto& [id, type] : nodes) {
        LOG_INFO(" Node {}: {}", id, (type == ENodeType::INPUT ? "INPUT" :
                                        type == ENodeType::OUTPUT ? "OUTPUT" :
                                        type == ENodeType::BIAS ? "BIAS" : "HIDDEN"));
    }

    // Now validate the complete network
    rebuildNetwork();
    validate();
}

void Genome::mutateActivations() {
    for (auto& gene : genes) {
        gene.activation.mutate(config.activationConfig);
    }
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
    
    LOG_DEBUG("Finding possible connections:");
    LOG_DEBUG("Current nodes:");
    for (const auto& [id, type] : nodes) {
        LOG_DEBUG(" Node {}: {}", id, (type == ENodeType::INPUT ? "INPUT" :
                      type == ENodeType::OUTPUT ? "OUTPUT" :
                      type == ENodeType::BIAS ? "BIAS" : "HIDDEN"));
    }
    
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
            LOG_TRACE("  Found possible connection: {} -> {}", fromId, toId);
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
    markDirty();

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
    markDirty();
    
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

void Genome::addConnection(NodeId from, NodeId to, double weight, bool validateAfter, core::EActivationType activation) {
    markDirty();
    
    // Ensure both nodes exist before adding connection
    if (from != -1) {  // Skip check for bias node
        auto fromType = nodes.find(from) != nodes.end() ? nodes[from] : ENodeType::HIDDEN;
        ensureNodeExists(from, fromType);
    }
    
    auto toType = nodes.find(to) != nodes.end() ? nodes[to] : ENodeType::HIDDEN;
    ensureNodeExists(to, toType);
    
    // Now add the connection
    Gene gene(from, to, weight, true, InnovationTracker::getNextInnovation());
    gene.activation = ActivationGene(activation);
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
    markDirty();
    
    validate();

    LOG_DEBUG("Starting connection mutation with {} inputs and {} outputs",
        std::count_if(nodes.begin(), nodes.end(),
                  [](const auto& pair) { return pair.second == ENodeType::INPUT; }),
        std::count_if(nodes.begin(), nodes.end(),
                  [](const auto& pair) { return pair.second == ENodeType::OUTPUT; }));
              
    // Make a deep copy backup using copy constructor
    Genome backup(*this);
    LOG_TRACE("Created backup - validating...");
    backup.validate();

    try {
        LOG_TRACE("Finding possible connections...");
        auto possibleConnections = findPossibleConnections();
        if (possibleConnections.empty()) {
            LOG_TRACE("No possible connections found");
            return false;
        }
        
        // Print current node structure before mutation
        LOG_TRACE("Current node structure before mutation:");
        for (const auto& [id, type] : nodes) {
            LOG_TRACE("  Node {}: {}", id, (type == ENodeType::INPUT ? "INPUT" :
                          type == ENodeType::OUTPUT ? "OUTPUT" :
                          type == ENodeType::BIAS ? "BIAS" : "HIDDEN"));
        }
        
        // Select random possible connection
        std::uniform_int_distribution<size_t> dist(0, possibleConnections.size() - 1);
        auto [fromId, toId] = possibleConnections[dist(rng)];
        
        // Generate random weight
        std::uniform_real_distribution<double> weightDist(-config.newWeightRange, config.newWeightRange);
        double weight = weightDist(rng);
        
        LOG_DEBUG("Adding connection from {} to {} with weight {}", fromId, toId, weight);
                
        // Add new connection with validation
        addConnection(fromId, toId, weight);
        
        validate();

        return true;
    } catch (const std::exception& e) {
        LOG_WARN("Connection mutation failed: {}", e.what());
        LOG_WARN("Restoring from backup...");
        
        // Verify backup before restore
        backup.validate();
        
        // Restore using assignment operator
        *this = std::move(backup);
        
        // Verify restoration
        validate();
        
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
    markDirty();
    
    LOG_DEBUG("Starting node mutation");
    LOG_DEBUG("Current genome state:");
    LOG_DEBUG("  Inputs: {}", std::count_if(nodes.begin(), nodes.end(),
        [](const auto& pair) { return pair.second == ENodeType::INPUT; }));
    LOG_DEBUG("  Outputs: {}", std::count_if(nodes.begin(), nodes.end(),
        [](const auto& pair) { return pair.second == ENodeType::OUTPUT; }));
    LOG_DEBUG("  Total nodes: {}", nodes.size());
    LOG_DEBUG("  Total genes: {}", genes.size());
    
    // Make a backup before mutation
    auto backup = *this;
    
    try {
        // Get enabled genes
        std::vector<std::reference_wrapper<Gene>> enabledGenes;
        for (auto& gene : genes) {
            if (gene.enabled) {
                enabledGenes.push_back(gene);
            }
        }
        
        if (enabledGenes.empty()) {
            LOG_DEBUG("No enabled genes available for node mutation");
            return false;
        }
        
        std::uniform_int_distribution<size_t> geneDist(0, enabledGenes.size() - 1);
        Gene& selectedGene = enabledGenes[geneDist(rng)];
        
        LOG_TRACE("Selected gene: {} -> {}", selectedGene.inputNode, selectedGene.outputNode);
        
        // Store old connection info
        NodeId fromNode = selectedGene.inputNode;
        NodeId toNode = selectedGene.outputNode;
        double oldWeight = selectedGene.weight;

        // Disable selected gene
        selectedGene.enabled = false;
        
        // Add new node
        NodeId newNodeId = getNextNodeId();
        LOG_TRACE("Creating new node: {}", newNodeId);

        addNode(newNodeId, ENodeType::HIDDEN, false);
        
        // Add new connections
        LOG_TRACE("Adding new connections");
        addConnection(fromNode, newNodeId, 1.0, false);  // Weight 1.0 to the new node
        addConnection(newNodeId, toNode, oldWeight, false);  // Keep old weight to output
        
        LOG_DEBUG("Node mutation complete. New structure:");
        for (const auto& [id, type] : nodes) {
            LOG_DEBUG("  Node {}: {}", id, (type == ENodeType::INPUT ? "INPUT" :
                          type == ENodeType::OUTPUT ? "OUTPUT" :
                          type == ENodeType::BIAS ? "BIAS" : "HIDDEN"));
        }

        // Validate the modified network
        validate();

        LOG_DEBUG("Node mutation successful");
        return true;
        
    } catch (const std::exception& e) {
        LOG_WARN("Node mutation failed: {}", e.what());
        // Restore backup on failure
        *this = backup;
        return false;
    }
}

void Genome::mutateWeights() {
    markDirty();
    
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

double Genome::compatibilityDistance(const Genome& genome1, const Genome& genome2) {
    const auto& genes1 = genome1.getGenes();
    const auto& genes2 = genome2.getGenes();
    
    double weightDiff = 0.0;
    double actDiff = 0.0;
    int32_t matchingGenes = 0;
    int32_t disjointGenes = 0;
    
    size_t i = 0, j = 0;
    while (i < genes1.size() && j < genes2.size()) {
        if (genes1[i].innovation == genes2[j].innovation) {
            weightDiff += std::abs(genes1[i].weight - genes2[j].weight);
            
            // Calculate activation difference using existing config
            double actCompatibility = Gene::getActivationCompatibility(genes1[i], genes2[j]);
            actDiff += actCompatibility * genome1.getConfig().activationConfig.compatibilityWeight;
            
            matchingGenes++;
            i++;
            j++;
        }
        else if (genes1[i].innovation < genes2[j].innovation) {
            disjointGenes++;
            i++;
        }
        else {
            disjointGenes++;
            j++;
        }
    }
    
    disjointGenes += (genes1.size() - i) + (genes2.size() - j);
    
    const double c1 = 1.0;  // Disjoint coefficient
    const double c2 = 1.0;  // Weight coefficient
    
    double distance = (c1 * disjointGenes) / std::max(genes1.size(), genes2.size());
    if (matchingGenes > 0) {
        distance += c2 * (weightDiff / matchingGenes);
        distance += actDiff / matchingGenes;
    }
    
    return distance;
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

std::vector<double> Genome::activate(const std::vector<double>& inputs) {
    // Validate network before any operation
    validate();

    if (networkDirty) {
        size_t newHash = calculateTopologyHash();
        if (newHash != topologyHash) {
            topologyHash = newHash;
            activationCache.clear();
        }
        networkDirty = false;
    }

    // Check cache
    size_t inputHash = calculateInputHash(inputs);
    auto cacheIt = activationCache.find(inputHash);
    if (cacheIt != activationCache.end()) {
        return cacheIt->second;
    }

    // Perform network activation
    auto outputs = network->activate(inputs);
    
    // Cache result
    activationCache[inputHash] = outputs;
    pruneCache();
    
    return outputs;
}

bool Genome::validate() const {
    std::vector<std::string> errors;
    
    // Debug output for node structure
    LOG_DEBUG("Validating genome...");
    LOG_DEBUG("Node structure:");
    for (const auto& [id, type] : nodes) {
        LOG_DEBUG("  Node {}: {}", id, (type == ENodeType::INPUT ? "INPUT" :
                      type == ENodeType::OUTPUT ? "OUTPUT" :
                      type == ENodeType::BIAS ? "BIAS" : "HIDDEN"));
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
    LOG_DEBUG("Gene structure:");
    for (const auto& gene : genes) {
        LOG_DEBUG("  Gene: {} -> {} (weight: {}, enabled {})", gene.inputNode, gene.outputNode,
            gene.weight, (gene.enabled ? "yes" : "no"));
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
    std::stringstream ss;
    for (const auto& id : reachable) {
        ss << id << " ";
    }
    LOG_DEBUG("Reachable nodes: {}", ss.str());
    
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
        LOG_ERROR("Validation failed with errors: {}", errorMsg);
        throw std::runtime_error(errorMsg);
    }
    
    LOG_DEBUG("Validation successful");
    return true;
}

void Genome::markDirty() {
    networkDirty = true;
    activationCache.clear();
    topologyHash = 0;
}

size_t Genome::calculateTopologyHash() const {
    size_t hash = 0;
    for (const auto& [id, type] : nodes) {
        hash ^= std::hash<int32_t>{}(id) + 0x9e3779b9;
    }
    for (const auto& gene : genes) {
        if (gene.enabled) {
            size_t geneHash = std::hash<int32_t>{}(gene.inputNode) ^
                            (std::hash<int32_t>{}(gene.outputNode) << 1) ^
                            (std::hash<double>{}(gene.weight) << 2);
            hash ^= geneHash + 0x9e3779b9;
        }
    }
    return hash;
}

size_t Genome::calculateInputHash(const std::vector<double>& inputs) const {
    size_t hash = topologyHash;
    for (const auto& input : inputs) {
        hash ^= std::hash<double>{}(input) + 0x9e3779b9;
    }
    return hash;
}

void Genome::pruneCache() {
    if (activationCache.size() > MAX_CACHE_SIZE) {
        size_t targetSize = MAX_CACHE_SIZE / 2;
        while (activationCache.size() > targetSize) {
            activationCache.erase(activationCache.begin());
        }
    }
}

}
}
