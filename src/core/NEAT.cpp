#include "NEAT.hpp"
#include "InnovationTracker.hpp"
#include "Species.hpp"
#include <algorithm>
#include <numeric>

namespace neat {
namespace core {

NEAT::NEAT(const Config& config)
    : config(config)
    , population(config.inputSize, config.outputSize, config.populationConfig) {}

void NEAT::evolve(int32_t generations, 
                 const std::function<double(const Genome&)>& fitnessFunction,
                 const std::function<void(int32_t, const Population&)>& generationCallback) {
    
    for (int32_t gen = 0; gen < generations; ++gen) {
        // Evaluate population
        population.evolve(fitnessFunction);
        
        // Update statistics
        updateStats();
        
        // Callback if provided
        if (generationCallback) {
            generationCallback(gen, population);
        }
    }
}

const Genome& NEAT::getBestGenome() const {
    return population.getBestGenome();
}

const Population& NEAT::getPopulation() const {
    return population;
}

NEAT::Stats NEAT::getStats() const {
    return currentStats;
}

void NEAT::updateStats() {
    const auto& species = population.getSpecies();
    
    // Find best fitness
    currentStats.bestFitness = population.getBestGenome().getFitness();
    
    // Calculate average fitness
    double totalFitness = 0.0;
    int32_t totalGenomes = 0;
    
    for (const auto& spec : species) {
        for (const auto& genomeIdx : spec.getMembers()) {
            totalFitness += population.getGenome(genomeIdx).getFitness();
            totalGenomes++;
        }
    }
    
    currentStats.averageFitness = totalGenomes > 0 ? totalFitness / totalGenomes : 0.0;
    currentStats.numSpecies = species.size();
    currentStats.generation++;
}

}
}
