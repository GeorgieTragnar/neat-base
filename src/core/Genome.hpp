// genome.hpp
#pragma once
#include <cstdint>
#include <vector>
#include <map>

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
    };
    
    Genome() = default;
    explicit Genome(const Config& config);
    
    // Core genome operations
    void addNode(NodeId id, ENodeType type);
    void addGene(const Gene& gene);
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
    Fitness fitness = 0.0;
    Fitness adjustedFitness = 0.0;
    int32_t speciesId = -1;
    int32_t maxNodeId = -1;
    Config config;
    
    // Helper methods
    bool hasCycle() const;
    bool isValidConnection(NodeId from, NodeId to) const;
    std::vector<NodeId> getTopologicalOrder() const;
    void sanitizeNetwork();
    static double activation(double x);
};

}
}
