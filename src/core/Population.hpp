// population.hpp
#pragma once
#include <cstdint>
#include <functional>
#include <random>
#include "core/Genome.hpp"

namespace neat {
namespace core {

class Genome;
class Species;

class Population {
public:
    struct Config {
        int32_t populationSize = 150;
        double survivalThreshold = 0.2;
        double compatibilityThreshold = 3.0;
        int32_t speciesTargetSize = 5;
        int32_t tournamentSize = 3;
        double elitismRate = 0.1;
        Genome::Config genomeConfig;
        
        Config() = default;
        Config(int32_t inputs, int32_t outputs)
            : genomeConfig(inputs, outputs) {}
    };
    
    Population(int32_t inputSize, int32_t outputSize, const Config& config);
    
    void evolve(const std::function<double(const Genome&)>& fitnessFunction);
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

    std::mt19937 rng;
};

}
}
