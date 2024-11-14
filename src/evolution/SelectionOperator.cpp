#include "SelectionOperator.hpp"
#include "core/Population.hpp"
#include "core/Genome.hpp"
#include "core/Species.hpp"
#include <algorithm>
#include <numeric>

namespace neat {
namespace evolution {

SelectionOperator::SelectionOperator(const Config& config)
    : config(config)
    , rng(std::random_device{}())
    , dist(0.0, 1.0)
{}

const core::Genome& SelectionOperator::selectParent(const core::Population& population) {
    // First select species based on fitness
    int32_t speciesIdx = selectSpecies(population.getSpecies());
    const auto& selectedSpecies = population.getSpecies()[speciesIdx];
    
    // Then select genome from that species
    return selectFromSpecies(selectedSpecies, population);
}

const core::Genome& SelectionOperator::tournamentSelect(const std::vector<core::Genome>& candidates) {
    if (candidates.empty()) {
        throw std::runtime_error("Cannot select from empty candidate pool");
    }
    
    std::uniform_int_distribution<size_t> indexDist(0, candidates.size() - 1);
    const core::Genome* best = nullptr;
    double bestFitness = std::numeric_limits<double>::lowest();
    
    for (int32_t i = 0; i < config.tournamentSize; ++i) {
        const auto& competitor = candidates[indexDist(rng)];
        double fitness = competitor.getFitness();
        
        if (fitness > bestFitness) {
            best = &competitor;
            bestFitness = fitness;
        }
    }
    
    return best ? *best : candidates[0];
}

std::vector<std::weak_ptr<core::Genome>> SelectionOperator::selectElites(const core::Population& population) {
    std::vector<std::weak_ptr<core::Genome>> elites;
    size_t eliteCount = static_cast<size_t>(population.size() * config.elitismRate);
    
    std::vector<size_t> indices(population.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    std::sort(indices.begin(), indices.end(),
        [&](size_t a, size_t b) {
            return population.getGenome(a).getFitness() > population.getGenome(b).getFitness();
        });
    
    for (size_t i = 0; i < eliteCount && i < indices.size(); ++i) {
        elites.push_back(std::make_shared<core::Genome>(population.getGenome(indices[i])));
    }
    
    return elites;
}

int32_t SelectionOperator::selectSpecies(const std::vector<core::Species>& species) {
    if (species.empty()) {
        throw std::runtime_error("No species available for selection");
    }
    
    // Calculate total adjusted fitness
    double totalFitness = 0.0;
    std::vector<double> fitnesses;
    fitnesses.reserve(species.size());
    
    for (const auto& s : species) {
        double speciesFitness = s.getBestFitness();
        fitnesses.push_back(speciesFitness);
        totalFitness += speciesFitness;
    }
    
    // Select species using roulette wheel selection
    double selected = dist(rng) * totalFitness;
    double sum = 0.0;
    
    for (size_t i = 0; i < species.size(); ++i) {
        sum += fitnesses[i];
        if (sum > selected) {
            return static_cast<int32_t>(i);
        }
    }
    
    return static_cast<int32_t>(species.size() - 1);
}

const core::Genome& SelectionOperator::selectFromSpecies(
    const core::Species& species, 
    const core::Population& population) {
        
    std::vector<core::Genome> candidates;
    candidates.reserve(species.getMembers().size());
    
    for (auto idx : species.getMembers()) {
        candidates.push_back(population.getGenome(idx));
    }
    
    return tournamentSelect(candidates);
}

}
}
