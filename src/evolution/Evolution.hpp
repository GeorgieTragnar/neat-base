// evolution.hpp
#pragma once

#include "MutationOperator.hpp"
#include "CrossoverOperator.hpp"
#include "SelectionOperator.hpp"
#include "SpeciesManager.hpp"
#include "FitnessEvaluator.hpp"
#include "core/Population.hpp"

namespace neat {
namespace evolution {

class Evolution {
public:
    struct Config {
        MutationOperator::Config mutationConfig;
        CrossoverOperator::Config crossoverConfig;
        SelectionOperator::Config selectionConfig;
        SpeciesManager::Config speciesConfig;
        FitnessEvaluator::Config fitnessConfig;
        
        int32_t populationSize = 150;
        int32_t maxGenerations = 100;
        double targetFitness = 1.0;
        bool stopOnTarget = true;
    };

    explicit Evolution(const Config& config);
    
    void initializePopulation(int32_t inputSize, int32_t outputSize);
    void evolve(const FitnessEvaluator::FitnessFunction& fitnessFunc);
    
    struct Stats {
        double bestFitness;
        double averageFitness;
        int32_t generation;
        int32_t speciesCount;
        const std::shared_ptr<core::Genome> bestGenome;
    };

    const Stats& getStats() const;
    const core::Population& getPopulation() const;

private:
    void epochComplete();
    bool shouldTerminate() const;
    void updateStats();

    Config config;
    core::Population population;
    Stats currentStats;
    
    MutationOperator mutationOp;
    CrossoverOperator crossoverOp;
    SelectionOperator selectionOp;
    SpeciesManager speciesManager;
    FitnessEvaluator fitnessEvaluator;
    
    std::mt19937_64 rng;
};

} // namespace evolution
} // namespace neat
