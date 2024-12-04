// genome.hpp
#pragma once

#include <random>
#include <set>
#include "network/Network.hpp"
#include "core/ActivationGene.hpp"
#include <iostream>

namespace neat {
namespace core {

enum class ENodeType;
class Gene;

class Genome {
public:
    using NodeId = int32_t;
    using Fitness = double;
    
    struct Config {
        // Core genome parameters
        int32_t inputSize;
        int32_t outputSize;
        double weightMutationRate;
        double weightPerturbationRate; 
        double newNodeRate;
        double newLinkRate;
        double weightPerturbationRange;
        double newWeightRange;
        int32_t maxHiddenNodes;

        ActivationGene::Config activationConfig;

        // Factory methods for child configs
        network::Network::Config createNetworkConfig() const {
            return network::Network::Config{
                inputSize,
                outputSize,
                activationConfig,
                false,  // allowRecurrent
                1.0     // biasValue
            };
        }

    public:
        Config() = delete;
        Config(int32_t inputs, int32_t outputs,
               double weightMutRate,
               double weightPerturbRate,
               double nodeRate,
               double linkRate,
               double perturbRange,
               double newWeightR,
               int32_t maxHidden,
               ActivationGene::Config actConfig)
            : inputSize(inputs)
            , outputSize(outputs)
            , weightMutationRate(weightMutRate)
            , weightPerturbationRate(weightPerturbRate)
            , newNodeRate(nodeRate)
            , newLinkRate(linkRate)
            , weightPerturbationRange(perturbRange)
            , newWeightRange(newWeightR)
            , maxHiddenNodes(maxHidden)
            , activationConfig(actConfig) {}
    };

    explicit Genome(const Config& config)
        : config(config)
        , network(std::make_unique<network::Network>(config.createNetworkConfig()))
        , fitness(0.0)
        , adjustedFitness(0.0)
        , speciesIdx(-1)
        , maxNodeIdx(-1)
        , networkDirty(true)
        , topologyHash(0)
        , rng(std::random_device{}()) {
        initMinimalTopology(config.inputSize, config.outputSize);
    }

    // Copy constructor
    Genome(const Genome& other)
        : genes(other.genes)
        , nodes(other.nodes)  // Copy nodes first
        , config(other.config)
        , fitness(other.fitness)
        , adjustedFitness(other.adjustedFitness)
        , speciesIdx(other.speciesIdx)
        , maxNodeIdx(other.maxNodeIdx)
        , networkDirty(true)
        , topologyHash(0)
        , network(std::make_unique<network::Network>(config.createNetworkConfig()))
        , rng(std::random_device{}()) {
        
        rebuildNetwork();  // Use existing method to rebuild network
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
        , networkDirty(other.networkDirty)
        , topologyHash(other.topologyHash)
        , config(other.config)
        , network(std::move(other.network))
        , rng(std::random_device{}()) {}
    
    // Copy assignment
    Genome& operator=(const Genome& other) {
        if (this != &other) {
            // Copy all members
            genes = other.genes;
            nodes = other.nodes;
            config = other.config;
            fitness = other.fitness;
            adjustedFitness = other.adjustedFitness;
            speciesIdx = other.speciesIdx;
            maxNodeIdx = other.maxNodeIdx;
            networkDirty = true;
            topologyHash = 0;
            activationCache.clear();
            
            rebuildNetwork();  // Use existing method to rebuild network
            
            validate();
        }
        return *this;
    }
    
    // Move assignment
    Genome& operator=(Genome&& other) noexcept {
        if (this != &other) {
            genes = std::move(other.genes);
            nodes = std::move(other.nodes);
            config = other.config;
            fitness = other.fitness;
            adjustedFitness = other.adjustedFitness;
            speciesIdx = other.speciesIdx;
            maxNodeIdx = other.maxNodeIdx;
            network = std::move(other.network);
        }
        return *this;
    }
    
    void setConfig(const Config& newConfig) { config = newConfig; }
    
    // Initialization
    void initMinimalTopology(int32_t inputs, int32_t outputs);
    NodeId getNextNodeId() noexcept { return ++maxNodeIdx; }

    // Core genome operations
    void addNode(NodeId id, ENodeType type, bool validateAfter = true);
    void addGene(const Gene& gene, bool validateAfter = true);
    void addConnection(NodeId from, NodeId to, double weight, bool validateAfter = true,
        core::EActivationType activation = core::EActivationType::SIGMOID);
    bool addConnectionMutation();
    bool addNodeMutation();
    void mutateWeights();
    void mutateActivations();
    
    // Network operations
    std::vector<double> activate(const std::vector<double>& inputs);
    bool validate() const;
    
    // Genetic operations
    static double compatibilityDistance(const Genome& genome1, const Genome& genome2);

    bool canAddGene(const Gene& gene) const;
    
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


    // Add these validation helper methods to the Genome class private section
    void ensureNodeExists(NodeId id, ENodeType type);
    std::vector<std::pair<NodeId, NodeId>> findPossibleConnections() const;
    bool wouldCreateCycle(NodeId fromId, NodeId toId) const;
    bool isValidGeneStructure(const Gene& gene) const;

    bool hasCycle() const;
    bool validateNodeStructure() const;
    bool validateGeneStructure() const;
    std::set<int32_t> getReachableNodes() const;
    void sanitizeNetwork();

    void markDirty();
    size_t calculateTopologyHash() const;
    size_t calculateInputHash(const std::vector<double>& inputs) const;
    void pruneCache();


    Config config;
    std::vector<Gene> genes;
    std::map<NodeId, ENodeType> nodes;
    std::unique_ptr<network::Network> network;
    Fitness fitness;
    Fitness adjustedFitness;
    int32_t speciesIdx;
    int32_t maxNodeIdx;
    bool networkDirty;
    size_t topologyHash;
    std::unordered_map<size_t, std::vector<double>> activationCache;
    std::mt19937 rng;

    static constexpr size_t MAX_CACHE_SIZE = 1000;
};

}
}
