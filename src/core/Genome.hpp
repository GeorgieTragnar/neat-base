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
    
    // Modifiers
    void setFitness(Fitness f) noexcept { fitness = f; }
    void setAdjustedFitness(Fitness f) noexcept { adjustedFitness = f; }
    void setSpecies(int32_t id) noexcept { speciesId = id; }

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
    
    void rebuildNetwork();
    bool isValidConnection(NodeId from, NodeId to) const;
};

}
}
