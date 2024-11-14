// genome.hpp
#pragma once

#include <random>
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
    
    Genome() = default;
    // Copy constructor
    Genome(const Genome& other)
        : genes(other.genes)
        , nodes(other.nodes)
        , fitness(other.fitness)
        , adjustedFitness(other.adjustedFitness)
        , speciesId(other.speciesId)
        , maxNodeId(other.maxNodeId)
        , config(other.config) {}
    
    // Move constructor
    Genome(Genome&& other) noexcept
        : genes(std::move(other.genes))
        , nodes(std::move(other.nodes))
        , fitness(other.fitness)
        , adjustedFitness(other.adjustedFitness)
        , speciesId(other.speciesId)
        , maxNodeId(other.maxNodeId)
        , config(other.config) {}
    
    // Copy assignment
    Genome& operator=(const Genome& other) {
        if (this != &other) {
            genes = other.genes;
            nodes = other.nodes;
            fitness = other.fitness;
            adjustedFitness = other.adjustedFitness;
            speciesId = other.speciesId;
            maxNodeId = other.maxNodeId;
            config = other.config;
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
            speciesId = other.speciesId;
            maxNodeId = other.maxNodeId;
            config = other.config;
        }
        return *this;
    }
    explicit Genome(const Config& config);
    
    // Initialization
    void initMinimalTopology(int32_t inputs, int32_t outputs);
    NodeId getNextNodeId() noexcept { return ++maxNodeId; }

    // Core genome operations
    void addNode(NodeId id, ENodeType type);
    void addGene(const Gene& gene);
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
    
    // Accessors
    Fitness getFitness() const noexcept { return fitness; }
    Fitness getAdjustedFitness() const noexcept { return adjustedFitness; }
    int32_t getSpecies() const noexcept { return speciesId; }
    const std::vector<Gene>& getGenes() const noexcept { return genes; }
    const std::map<NodeId, ENodeType>& getNodes() const noexcept { return nodes; }
    std::vector<Gene>& getGenes() noexcept { return genes; }
    std::map<NodeId, ENodeType>& getNodes() noexcept { return nodes; }
    
    // Modifiers
    void setFitness(Fitness f) noexcept { fitness = f; }
    void setAdjustedFitness(Fitness f) noexcept { adjustedFitness = f; }
    void setSpecies(int32_t id) noexcept { speciesId = id; }

    
    void rebuildNetwork();
    bool isValidConnection(NodeId from, NodeId to) const;

    Config getConfig() const { return config; }

private:
    std::vector<Gene> genes;
    std::map<NodeId, ENodeType> nodes;
    std::unique_ptr<network::Network> network;
    Fitness fitness = 0.0;
    Fitness adjustedFitness = 0.0;
    int32_t speciesId = -1;
    int32_t maxNodeId = -1;
    Config config;
    std::mt19937 rng;
};

}
}
