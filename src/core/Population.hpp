// population.hpp
#pragma once
#include <cstdint>
#include <functional>
#include <random>
#include "core/Genome.hpp"
#include "evolution/MutationOperator.hpp"

namespace neat {
namespace core {

class Genome;
class Species;

class Population {
public:
    struct Config {
        int32_t populationSize = 150;
        double survivalThreshold = 0.3;
        double compatibilityThreshold = 2.0;
        int32_t speciesTargetSize = 5;
        int32_t tournamentSize = 5;
        double elitismRate = 0.05;
        core::ActivationGene::Config activationConfig;
        Genome::Config genomeConfig;
        evolution::MutationOperator::Config mutationConfig;
        
        Config() = default;
        Config(int32_t inputs, int32_t outputs)
            : genomeConfig(inputs, outputs) {
                // Propagate activation config to mutation operator
                mutationConfig.activationConfig = activationConfig;
                genomeConfig.activationConfig = activationConfig;
                // TODO: need to redesign config propagation
            }
    };
    
    Population(int32_t inputSize, int32_t outputSize, const Config& config);
    
    void evolve(const std::function<double(Genome&)>& fitnessFunction);
    const Genome& getBestGenome() const;
    const std::vector<Species>& getSpecies() const;
    
    const Genome& getGenome(size_t idx) const { return genomes[idx]; }
    Genome& getGenome(size_t idx) { return genomes[idx]; }
    
    size_t size() const { return genomes.size(); }
    
private:
    std::vector<Genome> genomes;
    std::vector<Species> species;
    Config config;
    
    void speciate();
    void adjustFitness();
    Genome& selectParent();
    int32_t selectSpecies();
    void removeStaleSpecies();
    void removeWeakSpecies();

    evolution::MutationOperator mutationOp;
    std::mt19937 rng;
};

}
}
