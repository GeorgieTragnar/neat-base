#include "Evolution.hpp"
#include <algorithm>
#include <numeric>
#include <random>
#include "core/Population.hpp"

namespace neat {
namespace evolution {

Evolution::Evolution(const Config& config, core::Population& population)
    : config(config)
    , population(population)
    , mutationOp(config.mutationConfig)
    , crossoverOp(config.crossoverConfig)
    , selectionOp(config.selectionConfig)
    , speciesManager(config.speciesConfig)
    , fitnessEvaluator(config.fitnessConfig)
    , rng(std::random_device{}())
{
    currentStats = Stats{
        .bestFitness = 0.0,
        .averageFitness = 0.0,
        .generation = 0,
        .speciesCount = 0,
        .bestGenome = nullptr
    };
}

void Evolution::initializePopulation(int32_t inputSize, int32_t outputSize) {
    // Population already initialized in constructor
    speciesManager.speciate(population);
    updateStats();
}

void Evolution::evolve(const FitnessEvaluator::FitnessFunction& fitnessFunc) {
    while (!shouldTerminate()) {
        // Evaluate fitness
        fitnessEvaluator.evaluate(population, fitnessFunc);
        
        // Update species and adjust fitness
        speciesManager.updateSpecies(population);
        fitnessEvaluator.adjustFitness(population);
        
        // Create next generation
        std::vector<core::Genome> nextGen;
        nextGen.reserve(population.getConfig().populationSize);

        // Elitism - carry over best performers
        auto elites = selectionOp.selectElites(population);
        for (const auto& elite : elites) {
            if (auto genome = elite.lock()) {
                nextGen.emplace_back(*genome);
            }
        }

        // Fill rest of population through selection and reproduction
        while (nextGen.size() < population.getConfig().populationSize) {
            // Select parents
            const auto& parent1 = selectionOp.selectParent(population);
            const auto& parent2 = selectionOp.selectParent(population);
            
            // Perform crossover
            auto child = crossoverOp.crossover(parent1, parent2);
            
            // Mutate child
            mutationOp.mutate(child);
            
            nextGen.push_back(std::move(child));
        }

        // Safety check for population size
        if (nextGen.size() != population.size()) {
            throw std::runtime_error(
                "Population size mismatch: next generation size (" + 
                std::to_string(nextGen.size()) + 
                ") differs from current population size (" + 
                std::to_string(population.size()) + ")"
            );
        }

        // Replace old population one by one since there's no bulk replace
        for (size_t i = 0; i < population.size(); i++) {
            population.getGenome(i) = std::move(nextGen[i]);
        }
        
        // Remove stagnant species and adjust speciation
        speciesManager.removeStagnantSpecies(population);
        speciesManager.speciate(population);
        
        epochComplete();
    }
}

void Evolution::epochComplete() {
    currentStats.generation++;
    updateStats();
}

bool Evolution::shouldTerminate() const {
    if (currentStats.generation >= config.maxGenerations) {
        return true;
    }
    
    if (config.stopOnTarget && currentStats.bestFitness >= config.targetFitness) {
        return true;
    }
    
    return false;
}

void Evolution::updateStats() {
    double bestFitness = 0.0;
    double totalFitness = 0.0;
    const core::Genome* bestGenomePtr = nullptr;

    for (size_t i = 0; i < population.size(); i++) {
        const auto& genome = population.getGenome(i);
        double fitness = genome.getFitness();
        
        totalFitness += fitness;
        if (fitness > bestFitness) {
            bestFitness = fitness;
            bestGenomePtr = &genome;
        }
    }

    currentStats.bestFitness = bestFitness;
    currentStats.averageFitness = totalFitness / population.size();
    if (bestGenomePtr) {
        currentStats.bestGenome = std::make_shared<core::Genome>(*bestGenomePtr);
    }
    currentStats.speciesCount = speciesManager.getSpecies().size();
}

const Evolution::Stats& Evolution::getStats() const {
    return currentStats;
}

const core::Population& Evolution::getPopulation() const {
    return population;
}

}
}
