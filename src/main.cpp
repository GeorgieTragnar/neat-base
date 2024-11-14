#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <queue>
#include <unordered_map>
#include <functional>
#include <sstream>

// Random number generation
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(-1.0, 1.0);

enum class ENodeType {
    INPUT,
    HIDDEN,
    OUTPUT,
    BIAS
};

inline std::ostream& operator<<(std::ostream& os, const ENodeType& type) {
    switch (type) {
        case ENodeType::INPUT: return os << "input";
        case ENodeType::HIDDEN: return os << "hidden";
        case ENodeType::OUTPUT: return os << "output";
        case ENodeType::BIAS: return os << "bias";
        default: return os << "unknown";
    }
}

struct Gene {
    int inputNode;
    int outputNode;
    double weight;
    bool enabled;
    int innovation;
    
    Gene() : inputNode(0), outputNode(0), weight(0), enabled(true), innovation(0) {}
    
    Gene(int input, int output, double w, bool en, int innov)
        : inputNode(input), outputNode(output), weight(w), enabled(en), innovation(innov) {}
};

struct Species {
    std::vector<int> genomeIndices;
    double bestFitness;
    int staleness;
    Gene representative;
    double totalAdjustedFitness;  // Added this
    
    Species() : bestFitness(0), staleness(0), totalAdjustedFitness(0) {}
};

class Genome {
private:
    std::vector<Gene> genes;
    std::map<int, ENodeType> nodes; // node_id -> type ("input", "hidden", "output")
    double fitness;
    double adjustedFitness;
    int speciesIdx;
	int maxNodeId;
	static int nextNodeId;  // Add static counter for node IDs
    
    // Helper method to convert ENodeType to string (for debugging)
    static std::string ENodeTypeToString(ENodeType type) {
        switch (type) {
            case ENodeType::INPUT: return "input";
            case ENodeType::HIDDEN: return "hidden";
            case ENodeType::OUTPUT: return "output";
            case ENodeType::BIAS: return "bias";
            default: return "unknown";
        }
    }

    bool hasCycle() const {
        std::map<int, std::vector<int>> adjacencyList;
        std::set<int> allNodes;
        
        // Build adjacency list from enabled genes
        for (const auto& gene : genes) {
            if (gene.enabled) {
                adjacencyList[gene.inputNode].push_back(gene.outputNode);
                allNodes.insert(gene.inputNode);
                allNodes.insert(gene.outputNode);
            }
        }
        
        // For each node, check if it can reach itself
        for (int startNode : allNodes) {
            std::set<int> visited;
            std::set<int> recursionStack;
            
            std::function<bool(int)> hasCycleUtil = [&](int node) {
                visited.insert(node);
                recursionStack.insert(node);
                
                for (int neighbor : adjacencyList[node]) {
                    if (visited.find(neighbor) == visited.end()) {
                        if (hasCycleUtil(neighbor)) return true;
                    } else if (recursionStack.find(neighbor) != recursionStack.end()) {
                        return true;
                    }
                }
                
                recursionStack.erase(node);
                return false;
            };
            
            if (visited.find(startNode) == visited.end()) {
                if (hasCycleUtil(startNode)) return true;
            }
        }
        
        return false;
    }
    
    bool isValidConnection(int from, int to) const {
        // Check if nodes exist
        if (from != -1 && nodes.find(from) == nodes.end()) return false;
        if (nodes.find(to) == nodes.end()) return false;
        
        // Check node types
        if (from != -1 && nodes.at(from) == ENodeType::OUTPUT) return false;
        if (nodes.at(to) == ENodeType::INPUT) return false;
        
        // Check for existing connection
        for (const auto& gene : genes) {
            if (gene.enabled && gene.inputNode == from && gene.outputNode == to) {
                return false;
            }
        }
        
        return !wouldCreateCycle(from, to);
    }
    
    void sanitizeNetwork() {
        // First, disable any genes that reference non-existent nodes
        for (auto& gene : genes) {
            if (gene.enabled) {
                if (gene.inputNode != -1 && nodes.find(gene.inputNode) == nodes.end()) {
                    gene.enabled = false;
                    continue;
                }
                if (nodes.find(gene.outputNode) == nodes.end()) {
                    gene.enabled = false;
                    continue;
                }
            }
        }
        
        // Keep disabling genes until there are no cycles
        while (hasCycle()) {
            // Find a gene to disable (prefer newer genes)
            for (auto it = genes.rbegin(); it != genes.rend(); ++it) {
                if (it->enabled) {
                    it->enabled = false;
                    break;
                }
            }
        }
        
        // Ensure all nodes are reachable from inputs
        auto reachable = getReachableNodes();
        for (auto& gene : genes) {
            if (gene.enabled &&
                (reachable.find(gene.inputNode) == reachable.end() ||
                 reachable.find(gene.outputNode) == reachable.end())) {
                gene.enabled = false;
            }
        }
    }
    
	std::set<int> getReachableNodes() const {
		std::set<int> reachable;
		std::queue<int> queue;
		
		// Start with input nodes and bias
		for (const auto& [id, type] : nodes) {
			if (type == ENodeType::INPUT) {  // Changed from string comparison
				queue.push(id);
				reachable.insert(id);
			}
		}
		queue.push(-1); // bias node
		reachable.insert(-1);
		
		// Rest of the method remains the same
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

    // Add this method to validate the gene structure
    void validateGenes() const {
        for (const auto& gene : genes) {
            if (gene.enabled) {
                // Check input node
                if (gene.inputNode != -1 && nodes.find(gene.inputNode) == nodes.end()) {
                    std::stringstream ss;
                    ss << "Invalid input node in gene: " << gene.inputNode 
                       << " (max node ID: " << maxNodeId << ")";
                    throw std::runtime_error(ss.str());
                }
                
                // Check output node
                if (nodes.find(gene.outputNode) == nodes.end()) {
                    std::stringstream ss;
                    ss << "Invalid output node in gene: " << gene.outputNode 
                       << " (max node ID: " << maxNodeId << ")";
                    throw std::runtime_error(ss.str());
                }
            }
        }
    }
	
    std::vector<int> getTopologicalOrder() const {
		// Create a map of in-degrees
		std::map<int, int> inDegree;
		
		// Initialize in-degrees for all nodes
		for (const auto& [nodeId, ENodeType] : nodes) {
			inDegree[nodeId] = 0;
		}
		
		// Calculate in-degrees
		for (const auto& gene : genes) {
			if (gene.enabled) {
				if (inDegree.find(gene.outputNode) == inDegree.end()) {
					throw std::runtime_error("Invalid output node in gene: " + 
						std::to_string(gene.outputNode));
				}
				inDegree[gene.outputNode]++;
			}
		}
		
		// Initialize queue with nodes that have no incoming edges
		std::vector<int> queue;
		for (const auto& [nodeId, degree] : inDegree) {
			if (degree == 0) {
				queue.push_back(nodeId);
			}
		}
		
		if (queue.empty()) {
			throw std::runtime_error("No starting nodes found in topological sort");
		}
		
		// Process the queue
		std::vector<int> result;
		while (!queue.empty()) {
			int node = queue.front();
			queue.erase(queue.begin());
			result.push_back(node);
			
			// Reduce in-degrees of neighbors
			for (const auto& gene : genes) {
				if (gene.enabled && gene.inputNode == node) {
					inDegree[gene.outputNode]--;
					if (inDegree[gene.outputNode] == 0) {
						queue.push_back(gene.outputNode);
					}
				}
			}
		}
		
		// Check if we visited all nodes
		if (result.size() != nodes.size()) {
			std::cerr << "Warning: Topological sort didn't visit all nodes\n";
			std::cerr << "Visited " << result.size() << " nodes out of " << nodes.size() << "\n";
			
			// Print nodes that weren't visited
			std::set<int> visited(result.begin(), result.end());
			for (const auto& [nodeId, ENodeType] : nodes) {
				if (visited.find(nodeId) == visited.end()) {
					std::cerr << "Node " << nodeId << " (" << ENodeType << ") was not visited\n";
				}
			}
			
			throw std::runtime_error("Cycle detected or disconnected nodes in network");
		}
		
		return result;
	}

public:
    Genome() : fitness(0.0), adjustedFitness(0.0), speciesIdx(-1), maxNodeId(-1) {
        nodes[-1] = ENodeType::BIAS;
	}

	Genome(const Genome& other) 
		: genes(other.genes)
		, nodes(other.nodes)
		, fitness(other.fitness)
		, adjustedFitness(other.adjustedFitness)
		, speciesIdx(other.speciesIdx) 
		, maxNodeId(other.maxNodeId)
		{
			validateNetwork();
		}

	Genome& operator=(const Genome& other) {
		if (this != &other) {
			genes = other.genes;
			nodes = other.nodes;
			fitness = other.fitness;
			adjustedFitness = other.adjustedFitness;
			speciesIdx = other.speciesIdx;
			maxNodeId = other.maxNodeId;
			validateNetwork();
		}
		return *this;
	}

    // Modify validateNetwork to be more robust
    void validateNetwork() const {
        // First validate all nodes are properly mapped
        std::set<int> nodeIds;
        for (const auto& [id, type] : nodes) {
            if (id > maxNodeId && id != -1) {  // Allow bias node
                throw std::runtime_error("Node ID " + std::to_string(id) + 
                    " exceeds maxNodeId " + std::to_string(maxNodeId));
            }
            nodeIds.insert(id);
        }
        
        // Then validate all genes
        validateGenes();
        
        // Check for cycles
        if (hasCycle()) {
            throw std::runtime_error("Network contains cycles");
        }
    }

	// Add clone method for explicit copying
    Genome clone() const {
        Genome copy(*this);
        copy.validateNetwork();
        return copy;
    }

    void addGene(const Gene& gene) {
        genes.push_back(gene);
    }
    
    void addNode(int id, const ENodeType& type) {
		nodes[id] = type;
		if (id > maxNodeId && id != -1) {  // Don't count bias node
			maxNodeId = id;
		}
    }
    
    double activation(double x) const {
        // Steeper sigmoid for better XOR performance
        return 1.0 / (1.0 + std::exp(-4.9 * x));
    }

	int getMaxNodeId() const {
		return maxNodeId;
	}

    std::vector<double> feedForward(const std::vector<double>& inputs) const {
        try {
            // Validate inputs
            int inputCount = 0;
            int outputCount = 0;
            for (const auto& [nodeId, ENodeType] : nodes) {
            if (ENodeType == ENodeType::INPUT && nodeId != -1) inputCount++;  // Don't count bias node
            if (ENodeType == ENodeType::OUTPUT) outputCount++;
            }
            
            // Subtract 1 from inputCount to account for bias node
            if (inputs.size() != inputCount) {
                throw std::runtime_error("Invalid input size: got " + 
                    std::to_string(inputs.size()) + " expected " + 
                    std::to_string(inputCount - 1));
            }
            
            // Initialize nodeValues with input values and bias
            std::unordered_map<int, double> nodeValues;
            nodeValues[-1] = 1.0;  // Bias node
            
            int inputIdx = 0;
            for (const auto& [nodeId, ENodeType] : nodes) {
                if (ENodeType == ENodeType::INPUT && nodeId != -1) {
                    if (inputIdx >= inputs.size()) {
                        throw std::runtime_error("Input index out of bounds");
                    }
                    nodeValues[nodeId] = inputs[inputIdx++];
                }
            }
            
            // Process nodes in topological order
            auto sortedNodes = getTopologicalOrder();
            
            for (int node : sortedNodes) {
                auto nodeIt = nodes.find(node);
                if (nodeIt == nodes.end()) {
                    throw std::runtime_error("Node " + std::to_string(node) + " not found");
                }
                
                // Skip input nodes and bias node
                if (nodeIt->second == ENodeType::INPUT || node == -1) continue;
                
                double sum = 0.0;
                for (const auto& gene : genes) {
                    if (gene.enabled && gene.outputNode == node) {
                        auto it = nodeValues.find(gene.inputNode);
                        if (it == nodeValues.end()) {
                            throw std::runtime_error("No value found for node " + 
                                std::to_string(gene.inputNode));
                        }
                        sum += gene.weight * it->second;
                    }
                }
                
                nodeValues[node] = activation(sum);
            }
            
            // Collect outputs in order
            std::vector<double> outputs;
            for (const auto& [nodeId, ENodeType] : nodes) {
                if (ENodeType == ENodeType::OUTPUT) {
                    auto it = nodeValues.find(nodeId);
                    if (it == nodeValues.end()) {
                        throw std::runtime_error("No value computed for output node " + 
                            std::to_string(nodeId));
                    }
                    outputs.push_back(it->second);
                }
            }
            
            if (outputs.empty()) {
                throw std::runtime_error("No output values produced");
            }
            
            return outputs;
            
        } catch (const std::exception& e) {
            std::cerr << "Error in feedForward: " << e.what() << std::endl;
            throw;
        }
    }

    static double compatibilityDistance(const Genome& genome1, const Genome& genome2) {
        const double c1 = 1.0; // Excess weight
        const double c2 = 1.0; // Disjoint weight
        const double c3 = 0.4; // Weight diff weight
        
        std::map<int, const Gene*> genes1;
        std::map<int, const Gene*> genes2;
        
        for (const auto& gene : genome1.genes) genes1[gene.innovation] = &gene;
        for (const auto& gene : genome2.genes) genes2[gene.innovation] = &gene;
        
        int maxInnov1 = genes1.empty() ? 0 : genes1.rbegin()->first;
        int maxInnov2 = genes2.empty() ? 0 : genes2.rbegin()->first;
        
        int excess = 0;
        int disjoint = 0;
        double weightDiff = 0.0;
        int matching = 0;
        
        for (const auto& pair : genes1) {
            int innovation = pair.first;
            if (genes2.find(innovation) != genes2.end()) {
                matching++;
                weightDiff += std::abs(pair.second->weight - genes2[innovation]->weight);
            } else if (innovation > maxInnov2) {
                excess++;
            } else {
                disjoint++;
            }
        }
        
        for (const auto& pair : genes2) {
            int innovation = pair.first;
            if (genes1.find(innovation) == genes1.end() && innovation > maxInnov1) {
                excess++;
            } else if (genes1.find(innovation) == genes1.end()) {
                disjoint++;
            }
        }
        
        double n = std::max(genes1.size(), genes2.size());
        if (n < 20) n = 1;
        
        weightDiff = matching == 0 ? 0 : weightDiff / matching;
        return (c1 * excess + c2 * disjoint) / n + c3 * weightDiff;
    }
    
    void mutateWeights() {
        // More aggressive weight mutation
        for (auto& gene : genes) {
            if (dis(gen) < 0.9) { // Increased probability
                if (dis(gen) < 0.7) {
                    // Perturb weight
                    gene.weight += std::normal_distribution<>(0.0, 0.5)(gen);
                } else {
                    // Random new weight
                    gene.weight = std::normal_distribution<>(0.0, 2.0)(gen);
                }
            }
        }
    }

    bool wouldCreateCycle(int from, int to) const {
        // Create a temporary adjacency list
        std::map<int, std::vector<int>> adjacencyList;
        std::set<int> allNodes;
        
        // Add existing connections
        for (const auto& gene : genes) {
            if (gene.enabled) {
                adjacencyList[gene.inputNode].push_back(gene.outputNode);
                allNodes.insert(gene.inputNode);
                allNodes.insert(gene.outputNode);
            }
        }
        
        // Add the new connection
        adjacencyList[from].push_back(to);
        allNodes.insert(from);
        allNodes.insert(to);
        
        // Check for cycles
        std::set<int> visited;
        std::set<int> recursionStack;
        
        std::function<bool(int)> hasCycleUtil = [&](int node) {
            visited.insert(node);
            recursionStack.insert(node);
            
            for (int neighbor : adjacencyList[node]) {
                if (visited.find(neighbor) == visited.end()) {
                    if (hasCycleUtil(neighbor)) return true;
                } else if (recursionStack.find(neighbor) != recursionStack.end()) {
                    return true;
                }
            }
            
            recursionStack.erase(node);
            return false;
        };
        
        for (int node : allNodes) {
            if (visited.find(node) == visited.end()) {
                if (hasCycleUtil(node)) return true;
            }
        }
        
        return false;
    }

    bool addConnectionMutation(int& innovationNumber) {
        std::vector<std::pair<int, int>> possibleConnections;
        
        // Find all valid potential connections
        for (const auto& [id1, type1] : nodes) {
            if (type1 == ENodeType::OUTPUT) continue;
            
            for (const auto& [id2, type2] : nodes) {
                if (type2 == ENodeType::INPUT || id2 > 100) continue; // Limit node IDs
                if (id1 == id2) continue;
                
                if (isValidConnection(id1, id2)) {
                    possibleConnections.push_back({id1, id2});
                }
            }
        }
        
        if (possibleConnections.empty()) return false;
        
        // Choose random possible connection
        size_t idx = std::uniform_int_distribution<>(0, possibleConnections.size()-1)(gen);
        auto [input, output] = possibleConnections[idx];
        
        // Add new connection with limited innovation number
        innovationNumber = std::min(innovationNumber + 1, 1000);
        genes.push_back(Gene(input, output, 
            std::normal_distribution<>(0.0, 1.0)(gen), true, innovationNumber));
        
        return true;
    }

    // Add methods to access/modify node IDs
    static void resetNodeIds(int startId = 0) {
        nextNodeId = startId;
    }
    
    static int getNextNodeId() {
        return nextNodeId;
    }

    // Add a method to validate/fix the network
    void validate() {
        sanitizeNetwork();
    }

    bool addNodeMutation(int& innovationNumber) {
        if (genes.empty()) return false;

        // Get enabled genes
        std::vector<size_t> enabledGenes;
        for (size_t i = 0; i < genes.size(); i++) {
            if (genes[i].enabled) enabledGenes.push_back(i);
        }
        
        if (enabledGenes.empty()) return false;
        
        // Choose random enabled gene
        size_t geneIdx = enabledGenes[std::uniform_int_distribution<>(0, enabledGenes.size()-1)(gen)];
        Gene& chosen = genes[geneIdx];
        
        // Create new node ID
        int newNodeId = maxNodeId + 1;
        if (newNodeId > 100) return false; // Prevent runaway node creation
        
        // Verify that this node wouldn't create cycles
        if (chosen.inputNode != -1 && wouldCreateCycle(chosen.inputNode, newNodeId)) return false;
        if (wouldCreateCycle(newNodeId, chosen.outputNode)) return false;
        
        // Disable old connection
        chosen.enabled = false;
        
        // Add new node
        addNode(newNodeId, ENodeType::HIDDEN);
        maxNodeId = newNodeId;
        
        // Add new connections with limited innovation numbers
        innovationNumber = std::min(innovationNumber + 1, 1000); // Limit innovation numbers
        genes.push_back(Gene(chosen.inputNode, newNodeId, 1.0, true, innovationNumber));
        
        innovationNumber = std::min(innovationNumber + 1, 1000);
        genes.push_back(Gene(newNodeId, chosen.outputNode, chosen.weight, true, innovationNumber));
        
        return true;
    }

    void setFitness(double f) { fitness = f; }
    void setAdjustedFitness(double f) { adjustedFitness = f; }
    void setSpecies(int idx) { speciesIdx = idx; }
    
    double getFitness() const { return fitness; }
    double getAdjustedFitness() const { return adjustedFitness; }
    int getSpecies() const { return speciesIdx; }
    const std::vector<Gene>& getGenes() const { return genes; }
    std::vector<Gene>& getGenes() { return genes; }
    const std::map<int, ENodeType>& getNodes() const { return nodes; }  // Changed return type
};

int Genome::nextNodeId = 0;

class NEAT {
private:
    int inputSize;
    int outputSize;
    int populationSize;
    int innovationNumber;
	int globalMaxNodeId;
    double compatibilityThreshold = 4.0;
    std::vector<Genome> population;
    std::vector<Species> species;

    static constexpr double COMPATIBILITY_THRESHOLD = 3.0;
    static constexpr double SURVIVAL_THRESHOLD = 0.2;
    static constexpr int TOURNAMENT_SIZE = 3;
    
    
    // Add speciation parameters
    const double c1 = 1.0;  // Excess coefficient
    const double c2 = 1.0;  // Disjoint coefficient
    const double c3 = 0.4;  // Weight coefficient

	std::pair<double, double> calculateStats() {
        double bestFitness = -std::numeric_limits<double>::infinity();
        double totalFitness = 0;
        
        for (const auto& genome : population) {
            double fitness = genome.getFitness();
            bestFitness = std::max(bestFitness, fitness);
            totalFitness += fitness;
        }
        
        double avgFitness = totalFitness / population.size();
        return {bestFitness, avgFitness};
    }

	double compatibilityDistance(const Genome& genome1, const Genome& genome2) {
        const auto& genes1 = genome1.getGenes();
        const auto& genes2 = genome2.getGenes();
        
        std::map<int, const Gene*> innovations1;
        std::map<int, const Gene*> innovations2;
        
        for (const auto& gene : genes1) innovations1[gene.innovation] = &gene;
        for (const auto& gene : genes2) innovations2[gene.innovation] = &gene;
        
        int disjoint = 0;
        int excess = 0;
        double weightDiff = 0.0;
        int matching = 0;
        
        int maxInnov1 = genes1.empty() ? 0 : genes1.back().innovation;
        int maxInnov2 = genes2.empty() ? 0 : genes2.back().innovation;
        int maxInnov = std::max(maxInnov1, maxInnov2);
        
        // Count disjoint and excess genes
        for (const auto& gene : genes1) {
            if (innovations2.find(gene.innovation) == innovations2.end()) {
                if (gene.innovation > maxInnov2) excess++;
                else disjoint++;
            }
        }
        
        for (const auto& gene : genes2) {
            if (innovations1.find(gene.innovation) == innovations1.end()) {
                if (gene.innovation > maxInnov1) excess++;
                else disjoint++;
            }
        }
        
        // Calculate average weight difference of matching genes
        for (const auto& gene : genes1) {
            auto it = innovations2.find(gene.innovation);
            if (it != innovations2.end()) {
                weightDiff += std::abs(gene.weight - it->second->weight);
                matching++;
            }
        }
        
        double N = std::max(genes1.size(), genes2.size());
        if (N < 20) N = 1;
        
        weightDiff = matching == 0 ? 0 : weightDiff / matching;
        
        return (c1 * excess + c2 * disjoint) / N + c3 * weightDiff;
    }

    void speciate() {
        // Clear old species membership
        for (auto& s : species) {
            s.genomeIndices.clear();
        }
        
        // For each genome
        for (int i = 0; i < population.size(); i++) {
            bool found = false;
            
            // Try to find a species it's compatible with
            for (size_t j = 0; j < species.size(); j++) {
                if (!species[j].genomeIndices.empty()) {
                    // Compare against species representative
                    double dist = compatibilityDistance(population[i], 
                                                      population[species[j].genomeIndices[0]]);
                    
                    if (dist < compatibilityThreshold) {
                        species[j].genomeIndices.push_back(i);
                        population[i].setSpecies(j);
                        found = true;
                        break;
                    }
                }
            }
            
            // If no compatible species found, create a new one
            if (!found) {
                Species newSpecies;
                newSpecies.genomeIndices.push_back(i);
                population[i].setSpecies(species.size());
                species.push_back(newSpecies);
            }
        }
        
        // Remove empty species
        species.erase(
            std::remove_if(species.begin(), species.end(),
                [](const Species& s) { return s.genomeIndices.empty(); }),
            species.end()
        );
    }
    
    void adjustFitness() {
        for (auto& s : species) {
            for (int idx : s.genomeIndices) {
                population[idx].setAdjustedFitness(
                    population[idx].getFitness() / s.genomeIndices.size()
                );
            }
        }
    }
    
    Genome& selectParent() {
        // Sum all adjusted fitnesses
        double totalFitness = 0;
        for (const auto& genome : population) {
            totalFitness += genome.getAdjustedFitness();
        }
        
        // Select parent based on adjusted fitness
        double roll = dis(gen) * totalFitness;
        double sum = 0;
        for (auto& genome : population) {
            sum += genome.getAdjustedFitness();
            if (sum > roll) {
                return genome;
            }
        }
        return population.back();
    }
	
    int selectSpecies() {
        double totalFitness = 0;
        for (const auto& species : species) {
            totalFitness += species.totalAdjustedFitness;
        }
        
        double roll = dis(gen) * totalFitness;
        double sum = 0;
        
        for (size_t i = 0; i < species.size(); i++) {
            sum += species[i].totalAdjustedFitness;
            if (sum > roll) {
                return i;
            }
        }
        
        return species.size() - 1;
    }

    Genome selectParentFromSpecies(int speciesIdx) {
        const auto& currentSpecies = species[speciesIdx];
        
        // Tournament selection within species
        int bestIdx = -1;
        double bestFitness = -std::numeric_limits<double>::infinity();
        
        for (int i = 0; i < TOURNAMENT_SIZE && i < currentSpecies.genomeIndices.size(); i++) {
            int idx = currentSpecies.genomeIndices[
                std::uniform_int_distribution<>(0, currentSpecies.genomeIndices.size() - 1)(gen)
            ];
            
            double fitness = population[idx].getAdjustedFitness();
            if (fitness > bestFitness) {
                bestFitness = fitness;
                bestIdx = idx;
            }
        }
        
        // If tournament selection failed (shouldn't happen), return random member
        if (bestIdx == -1) {
            bestIdx = currentSpecies.genomeIndices[
                std::uniform_int_distribution<>(0, currentSpecies.genomeIndices.size() - 1)(gen)
            ];
        }
        
        return population[bestIdx];
    }
    
    // Modify crossover to use clone and validate
    Genome crossover(const Genome& parent1, const Genome& parent2) {
        try {
            Genome child;
            const Genome& primaryParent = parent1.getFitness() >= parent2.getFitness() ? parent1 : parent2;
            const Genome& secondaryParent = parent1.getFitness() >= parent2.getFitness() ? parent2 : parent1;
            
            // First validate parent node IDs
            int maxNodeId = std::max(primaryParent.getMaxNodeId(), secondaryParent.getMaxNodeId());
            
            // Copy nodes from primary parent and ensure node IDs are valid
            for (const auto& [nodeId, nodeType] : primaryParent.getNodes()) {
                if (nodeId != -1 && nodeId > maxNodeId) {
                    throw std::runtime_error("Invalid node ID in primary parent: " + std::to_string(nodeId));
                }
                child.addNode(nodeId, nodeType);
            }
            
            // Add matching genes
            std::map<int, const Gene*> primaryGenes, secondaryGenes;
            for (const auto& gene : primaryParent.getGenes()) primaryGenes[gene.innovation] = &gene;
            for (const auto& gene : secondaryParent.getGenes()) secondaryGenes[gene.innovation] = &gene;
            
            // First, add matching genes
            for (const auto& [innov, gene1] : primaryGenes) {
                auto it = secondaryGenes.find(innov);
                if (it != secondaryGenes.end()) {
                    // For matching genes, randomly choose one
                    const Gene* selectedGene = (dis(gen) < 0.5) ? gene1 : it->second;
                    
                    // Validate node IDs before adding
                    if (selectedGene->inputNode != -1 && selectedGene->inputNode > maxNodeId) continue;
                    if (selectedGene->outputNode > maxNodeId) continue;
                    
                    // Verify nodes exist in child
                    if (selectedGene->inputNode != -1 && 
                        child.getNodes().find(selectedGene->inputNode) == child.getNodes().end()) continue;
                    if (child.getNodes().find(selectedGene->outputNode) == child.getNodes().end()) continue;
                    
                    child.addGene(*selectedGene);
                } else {
                    // For disjoint/excess genes from primary parent
                    if (gene1->inputNode != -1 && gene1->inputNode > maxNodeId) continue;
                    if (gene1->outputNode > maxNodeId) continue;
                    
                    if (gene1->inputNode != -1 && 
                        child.getNodes().find(gene1->inputNode) == child.getNodes().end()) continue;
                    if (child.getNodes().find(gene1->outputNode) == child.getNodes().end()) continue;
                    
                    child.addGene(*gene1);
                }
            }
            
            // Check if child is valid
            bool hasEnabledGenes = false;
            for (const auto& gene : child.getGenes()) {
                if (gene.enabled) {
                    hasEnabledGenes = true;
                    break;
                }
            }
            
            if (!hasEnabledGenes) {
                return primaryParent.clone();
            }
            
            child.validateNetwork();
            return child;
            
        } catch (const std::exception& e) {
            std::cerr << "Error in crossover: " << e.what() << std::endl;
            return parent1.getFitness() >= parent2.getFitness() ? parent1.clone() : parent2.clone();
        }
    }
    
    void mutate(Genome& genome) {
        // Increase mutation rates
        if (dis(gen) < 0.8) genome.mutateWeights();
        if (dis(gen) < 0.1) genome.addNodeMutation(innovationNumber);  // Increased from 0.05
        if (dis(gen) < 0.15) genome.addConnectionMutation(innovationNumber);  // Increased from 0.05
    }

public:
    NEAT(int inputs, int outputs, int popSize)
        : inputSize(inputs)
        , outputSize(outputs)
        , populationSize(popSize)
        , innovationNumber(0)
        , compatibilityThreshold(3.0)
		, globalMaxNodeId(inputs + outputs - 1)
    {
        initializePopulation();
    }
    
    void initializePopulation() {
        for (int i = 0; i < populationSize; i++) {
            Genome genome;
            
            // Add input nodes
            for (int j = 0; j < inputSize; j++) {
                genome.addNode(j, ENodeType::INPUT);  // Changed from string to enum
            }
            
            // Add output nodes
            for (int j = 0; j < outputSize; j++) {
                genome.addNode(inputSize + j, ENodeType::OUTPUT);  // Changed from string to enum
            }
            
            // Add initial connections
            for (int in = -1; in < inputSize; in++) {
                for (int out = 0; out < outputSize; out++) {
                    innovationNumber++;
                    genome.addGene(Gene(in, inputSize + out, 
                        std::normal_distribution<>(0.0, 2.0)(gen), true, innovationNumber));
                }
            }
            
            population.push_back(genome);
        }
        
        // Initialize first species
        Species initialSpecies;
        for (int i = 0; i < populationSize; i++) {
            initialSpecies.genomeIndices.push_back(i);
            population[i].setSpecies(0);
        }
        species.push_back(initialSpecies);
    }

    void evolve(int generationCount, std::function<double(const Genome&)> fitnessFunc) {
        try {
            globalMaxNodeId = inputSize + outputSize - 1;  // Initialize global max node ID

            for (int generation = 0; generation < generationCount; generation++) {
                // Track max node ID across population
                for (const auto& genome : population) {
                    globalMaxNodeId = std::max(globalMaxNodeId, genome.getMaxNodeId());
                }
                
                // Evaluate fitness
                for (auto& genome : population) {
                    genome.validateNetwork();  // Validate before evaluation
                    double fitness = fitnessFunc(genome);
                    genome.setFitness(fitness);
                }
                
                // Dynamically adjust compatibility threshold to target ~10 species
                if (generation > 0 && generation % 10 == 0) {
                    if (species.size() > 15) compatibilityThreshold *= 1.2;
                    else if (species.size() < 5) compatibilityThreshold *= 0.8;
                }
                
                // Speciate the population
                speciate();
            
                // Remove stagnant species
                for (auto& s : species) {
                    double maxFitness = -std::numeric_limits<double>::infinity();
                    for (int idx : s.genomeIndices) {
                        maxFitness = std::max(maxFitness, population[idx].getFitness());
                    }
                    
                    if (maxFitness > s.bestFitness) {
                        s.bestFitness = maxFitness;
                        s.staleness = 0;
                    } else {
                        s.staleness++;
                    }
                }
                
                species.erase(
                    std::remove_if(species.begin(), species.end(),
                        [](const Species& s) { 
                            return s.staleness > 15 && s.genomeIndices.size() < 5; 
                        }),
                    species.end()
                );
                
                if (species.empty()) {
                    // If all species were removed, keep the best performers
                    species.push_back(Species());
                    std::sort(population.begin(), population.end(),
                        [](const Genome& a, const Genome& b) {
                            return a.getFitness() > b.getFitness();
                        });
                    for (int i = 0; i < std::min(5, (int)population.size()); i++) {
                        species[0].genomeIndices.push_back(i);
                    }
                }

                // Calculate adjusted fitness
                adjustFitness();
                
                // Calculate total adjusted fitness for each species
                for (auto& s : species) {
                    s.totalAdjustedFitness = 0;
                    for (int idx : s.genomeIndices) {
                        s.totalAdjustedFitness += population[idx].getAdjustedFitness();
                    }
                }
                
                // Create new population
                std::vector<Genome> newPopulation;
                newPopulation.reserve(populationSize);
                
                // Sort species by best fitness
                std::sort(species.begin(), species.end(),
                    [&](const Species& a, const Species& b) {
                        return a.bestFitness > b.bestFitness;
                    });
                    
                // Elitism: Keep best genome from each species
                for (const auto& s : species) {
                    if (s.genomeIndices.size() < 5) continue;
                    
                    int bestIdx = *std::max_element(s.genomeIndices.begin(), s.genomeIndices.end(),
                        [this](int a, int b) {
                            return population[a].getFitness() < population[b].getFitness();
                        });
                    
                    auto elite = population[bestIdx].clone();
                    elite.validateNetwork();
                    newPopulation.push_back(std::move(elite));
                }
                
                // Fill rest of new population
                while (newPopulation.size() < populationSize) {
                    int speciesIdx = selectSpecies();
                    
                    // If species is empty or index invalid, use first species
                    if (speciesIdx >= species.size() || species[speciesIdx].genomeIndices.empty()) {
                        speciesIdx = 0;
                    }

                    try {
                        // If still no valid species found, clone best genome
                        if (species[speciesIdx].genomeIndices.empty()) {
                            auto bestGenome = std::max_element(population.begin(), population.end(),
                                [](const Genome& a, const Genome& b) {
                                    return a.getFitness() < b.getFitness();
                                });
                            newPopulation.push_back(bestGenome->clone());
                            continue;
                        }
                        
                        Genome child;
                        
                        // Sexual reproduction
                        if (dis(gen) < 0.9 && species[speciesIdx].genomeIndices.size() > 1) {
                            Genome parent1 = selectParentFromSpecies(speciesIdx);
                            Genome parent2 = selectParentFromSpecies(speciesIdx);
                            child = crossover(parent1, parent2);
                        } else {
                            // Asexual reproduction
                            child = selectParentFromSpecies(speciesIdx).clone();
                        }
                        
                        // Mutations
                        if (dis(gen) < 0.8) child.mutateWeights();
                        
                        // Only add nodes if network isn't too complex
                        if (dis(gen) < 0.1 && child.getMaxNodeId() < 15) {
                            if (!child.addNodeMutation(innovationNumber)) {
                                // If node mutation fails, try connection mutation instead
                                child.addConnectionMutation(innovationNumber);
                            }
                        } else if (dis(gen) < 0.15) {
                            child.addConnectionMutation(innovationNumber);
                        }
                        
                        // Validate and add to population
                        try {
                            child.validateNetwork();
                            newPopulation.push_back(std::move(child));
                        } catch (const std::exception& e) {
                            // If validation fails, clone a parent
                            child = selectParentFromSpecies(speciesIdx).clone();
                            child.validateNetwork();
                            newPopulation.push_back(std::move(child));
                        }
                        
                    } catch (const std::exception& e) {
                        std::cerr << "Error creating child: " << e.what() << std::endl;
                        // On error, clone the best genome from the species
                        if (!species[speciesIdx].genomeIndices.empty()) {
                            auto bestIdx = *std::max_element(
                                species[speciesIdx].genomeIndices.begin(),
                                species[speciesIdx].genomeIndices.end(),
                                [this](int a, int b) {
                                    return population[a].getFitness() < population[b].getFitness();
                                });
                            newPopulation.push_back(population[bestIdx].clone());
                        }
                    }
                }
                
                // Validate new population before assignment
                for (const auto& genome : newPopulation) {
                    genome.validateNetwork();
                }
                
                population = std::move(newPopulation);
                
                // Optional: Print statistics
                auto [bestFitness, avgFitness] = calculateStats();
                if (generation % 10 == 0) {
                    std::cout << "Generation " << generation 
                            << ": Best Fitness = " << bestFitness 
                            << ", Average Fitness = " << avgFitness 
                            << ", Species = " << species.size() << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in evolution: " << e.what() << std::endl;
            throw;
        }
    }

    const Genome& getBestGenome() const {
        return population[0];
    }
};

int main() {
    NEAT neat(2, 1, 150);

    auto fitnessFunc = [](const Genome& genome) {
        std::vector<std::vector<double>> xorInputs = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
        std::vector<double> xorOutputs = {0, 1, 1, 0};
        
        double error = 0.0;
        for (size_t i = 0; i < xorInputs.size(); i++) {
            double output = genome.feedForward(xorInputs[i])[0];
            error += std::pow(output - xorOutputs[i], 2);
        }
        
        return 4.0 / (error + 1.0);  // Normalized fitness
    };

    neat.evolve(500, fitnessFunc);

    // Test best solution
    const Genome& bestGenome = neat.getBestGenome();
    std::cout << "\nTesting best solution:" << std::endl;
    std::vector<std::vector<double>> testInputs = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    for (const auto& inputs : testInputs) {
        double output = bestGenome.feedForward(inputs)[0];
        std::cout << "Input: [" << inputs[0] << ", " << inputs[1] 
                  << "], Output: " << output << std::endl;
    }

    return 0;
}