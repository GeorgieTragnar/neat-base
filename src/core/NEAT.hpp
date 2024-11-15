// neat.hpp
#pragma once
#include <cstdint>
#include "Population.hpp"
#include "Genome.hpp"

namespace neat {
namespace core {

class NEAT {
public:
    struct Config {
        int32_t inputSize;
        int32_t outputSize;
        Population::Config populationConfig;

        // Add explicit Genome config initialization
        Config(int32_t inputs, int32_t outputs) 
            : inputSize(inputs)
            , outputSize(outputs)
            , populationConfig(inputs, outputs) {}
    };
    
    explicit NEAT(const Config& config);
    
    void evolve(int32_t generations, 
                const std::function<double(const Genome&)>& fitnessFunction,
                const std::function<void(int32_t, const Population&)>& generationCallback = nullptr);
    
    const Genome& getBestGenome() const;
    const Population& getPopulation() const;
    
    // Statistics and monitoring
    struct Stats {
        double bestFitness;
        double averageFitness;
        int32_t numSpecies;
        int32_t generation;
    };
    
    Stats getStats() const;

private:
    Population population;
    Config config;
    Stats currentStats;
    
    void updateStats();
};

} // namespace core
} // namespace neat
