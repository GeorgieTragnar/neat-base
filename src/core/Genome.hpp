// genome.hpp
#pragma once

#include <random>
#include <set>
#include "network/Network.hpp"

namespace neat {
namespace core {

enum class ENodeType;
class Gene;

class Genome {
public:
    using NodeId = int32_t;
    using Fitness = double;
    
    struct Config {
        double weightMutationRate = 0.8;
        double weightPerturbationRate = 0.9;
        double newNodeRate = 0.03;
        double newLinkRate = 0.05;
        double weightPerturbationRange = 0.2;
        double newWeightRange = 2.0;
        network::Network::Config networkConfig;
    };
    
    // Enhanced constructors with validation
    Genome() 
        : config()  // Initialize with default config
        , network(std::make_unique<network::Network>(config.networkConfig))
        , fitness(0.0)
        , adjustedFitness(0.0)
        , speciesIdx(-1)
        , maxNodeIdx(-1)
        , rng(std::random_device{}()) {
        nodes[-1] = ENodeType::BIAS;
        network->addNode(-1, ENodeType::BIAS);  // Add bias to network
    }

    // Copy constructor
    Genome(const Genome& other)
        : genes(other.genes)
        , nodes(other.nodes)
        , fitness(other.fitness)
        , adjustedFitness(other.adjustedFitness)
        , speciesIdx(other.speciesIdx)
        , maxNodeIdx(other.maxNodeIdx)
        , config(other.config)
        , network(std::make_unique<network::Network>(config.networkConfig))
        , rng(std::random_device{}()) {
        rebuildNetwork();
        validate();
    }
    
    // Move constructor
    Genome(Genome&& other) noexcept
        : genes(std::move(other.genes))
        , nodes(std::move(other.nodes))
        , fitness(other.fitness)
        , adjustedFitness(other.adjustedFitness)
        , speciesIdx(other.speciesIdx)
        , maxNodeIdx(other.maxNodeIdx)
        , config(other.config)
        , network(std::move(other.network))
        , rng(std::random_device{}()) {}
    
    // Copy assignment
    Genome& operator=(const Genome& other) {
        if (this != &other) {
            genes = other.genes;
            nodes = other.nodes;
            fitness = other.fitness;
            adjustedFitness = other.adjustedFitness;
            speciesIdx = other.speciesIdx;
            maxNodeIdx = other.maxNodeIdx;
            config = other.config;
            network = std::make_unique<network::Network>(config.networkConfig);
            rebuildNetwork();
            validate();
        }
        return *this;
    }
    
    // Move assignment
    Genome& operator=(Genome&& other) noexcept {
        if (this != &other) {
            genes = std::move(other.genes);
            nodes = std::move(other.nodes);
            fitness = other.fitness;
            adjustedFitness = other.adjustedFitness;
            speciesIdx = other.speciesIdx;
            maxNodeIdx = other.maxNodeIdx;
            config = other.config;
            network = std::move(other.network);
        }
        return *this;
    }

    explicit Genome(const Config& config)
        : config(config)
        , fitness(0.0)
        , adjustedFitness(0.0)
        , speciesIdx(-1)
        , maxNodeIdx(-1)
        , network(std::make_unique<network::Network>(config.networkConfig))
        , rng(std::random_device{}()) {
            nodes[-1] = ENodeType::BIAS;  // Ensure bias node exists
            network->addNode(-1, ENodeType::BIAS);  // Add bias to network
        }
    
    // Initialization
    void initMinimalTopology(int32_t inputs, int32_t outputs);
    NodeId getNextNodeId() noexcept { return ++maxNodeIdx; }

    // Core genome operations
    void addNode(NodeId id, ENodeType type, bool validateAfter = true);
    void addGene(const Gene& gene, bool validateAfter = true);
    void addConnection(NodeId from, NodeId to, double weight);
    bool addConnectionMutation();
    bool addNodeMutation();
    void mutateWeights();
    
    // Network operations
    std::vector<double> activate(const std::vector<double>& inputs) const;
    bool validate() const;
    
    // Genetic operations
    static Genome crossover(const Genome& parent1, const Genome& parent2, const Config& config);
    static double compatibilityDistance(const Genome& genome1, const Genome& genome2);
    
    // Clone method for safe copying
    Genome clone() const;

    // Accessors
    Fitness getFitness() const noexcept { return fitness; }
    Fitness getAdjustedFitness() const noexcept { return adjustedFitness; }
    int32_t getSpecies() const noexcept { return speciesIdx; }
    const std::vector<Gene>& getGenes() const noexcept { return genes; }
    const std::map<NodeId, ENodeType>& getNodes() const noexcept { return nodes; }
    std::vector<Gene>& getGenes() noexcept { return genes; }
    std::map<NodeId, ENodeType>& getNodes() noexcept { return nodes; }
    
    // Modifiers
    void setFitness(Fitness f) noexcept { fitness = f; }
    void setAdjustedFitness(Fitness f) noexcept { adjustedFitness = f; }
    void setSpecies(int32_t id) noexcept { speciesIdx = id; }

    // Safe accessors and modifiers
    int getMaxNodeId() const { return maxNodeIdx; }
    void setMaxNodeId(int id) { maxNodeIdx = id; }
    
    void rebuildNetwork();
    bool isValidConnection(NodeId from, NodeId to) const;

    Config getConfig() const { return config; }

private:

    bool hasCycle() const;
    bool validateNodeStructure() const;
    bool validateGeneStructure() const;
    std::set<int32_t> getReachableNodes() const;
    void sanitizeNetwork();

    std::vector<Gene> genes;
    std::map<NodeId, ENodeType> nodes;
    std::unique_ptr<network::Network> network;
    Fitness fitness = 0.0;
    Fitness adjustedFitness = 0.0;
    int32_t speciesIdx = -1;
    int32_t maxNodeIdx = -1;
    static int32_t nextNodeIdx;
    Config config;
    std::mt19937 rng;
};

}
}
