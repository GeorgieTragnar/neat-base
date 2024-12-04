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
        
        // Remove populationSize since it's defined by NEAT's Population
        int32_t maxGenerations;
        double targetFitness;
        bool stopOnTarget;

        Config() = delete;
        Config(const MutationOperator::Config& mutConfig,
            const CrossoverOperator::Config& crossConfig,
            const SelectionOperator::Config& selectConfig,
            const SpeciesManager::Config& speciesConfig,
            const FitnessEvaluator::Config& fitnessConfig,
            int32_t maxGen,
            double targetFit,
            bool stopTarget)
            : mutationConfig(mutConfig)
            , crossoverConfig(crossConfig)
            , selectionConfig(selectConfig)
            , speciesConfig(speciesConfig)
            , fitnessConfig(fitnessConfig)
            , maxGenerations(maxGen)
            , targetFitness(targetFit)
            , stopOnTarget(stopTarget) {}
    };

    explicit Evolution(const Config& config, core::Population& population);
    
    void initializePopulation(int32_t inputSize, int32_t outputSize);
    void evolve(const FitnessEvaluator::FitnessFunction& fitnessFunc);
    
    struct Stats {
        double bestFitness;
        double averageFitness;
        int32_t generation;
        int32_t speciesCount;
        std::shared_ptr<core::Genome> bestGenome;
    };

    const Stats& getStats() const;
    const core::Population& getPopulation() const;

private:
    void epochComplete();
    bool shouldTerminate() const;
    void updateStats();

    Config config;
    core::Population& population;
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
