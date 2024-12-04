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
        int32_t maxGenerations;
        double targetFitness;
        bool stopOnTarget;

        // Factory method for Population config
        Population::Config createPopulationConfig() const {
            return Population::Config{
                inputSize,
                outputSize,
                500,  // Population size 
                0.2,  // Survival threshold
                3.0,  // Compatibility threshold
                10,   // Species target
                3,    // Tournament size
                0.1   // Elitism rate
            };
        }

        Config() = delete;
        Config(int32_t inputs, int32_t outputs)
            : inputSize(inputs)
            , outputSize(outputs)
            , maxGenerations(100)
            , targetFitness(1.0)
            , stopOnTarget(true) {}
    };
    
    explicit NEAT(const Config& config);
    
    void evolve(int32_t generations, 
                const std::function<double(Genome&)>& fitnessFunction,
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
    Config config;
    Population population;
    Stats currentStats;
    
    void updateStats();
};

} // namespace core
} // namespace neat
