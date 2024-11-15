#include "Population.hpp"
#include "Genome.hpp"
#include "Species.hpp"
#include <algorithm>
#include <random>
#include <iostream>

namespace neat {
namespace core {

Population::Population(int32_t inputSize, int32_t outputSize, const Config& config)
    : config(config)
    , rng(std::random_device{}()) {
    
    std::cout << "Creating population with " << inputSize << " inputs and " 
              << outputSize << " outputs" << std::endl;
              
    genomes.reserve(config.populationSize);
    
    // Create initial genomes with proper config
    for (int32_t i = 0; i < config.populationSize; ++i) {
        // Pass the genome config from population config
        Genome genome(config.genomeConfig);
        genome.initMinimalTopology(inputSize, outputSize);
        
        // Verify genome structure
        auto& nodes = genome.getNodes();
        int inputCount = 0, outputCount = 0;
        for (const auto& [id, type] : nodes) {
            if (type == ENodeType::INPUT) inputCount++;
            if (type == ENodeType::OUTPUT) outputCount++;
        }
        
        std::cout << "Genome " << i << " initialized with " 
                  << inputCount << " inputs and "
                  << outputCount << " outputs" << std::endl;
                  
        genomes.push_back(std::move(genome));
    }
    
    std::cout << "Population initialized with " << genomes.size() 
              << " genomes" << std::endl;
              
    speciate();
}

void Population::evolve(const std::function<double(const Genome&)>& fitnessFunction) {
    // Evaluate fitness
    for (auto& genome : genomes) {
        genome.setFitness(fitnessFunction(genome));
    }
    
    // Update species and adjust fitness
    speciate();
    adjustFitness();
    removeStaleSpecies();
    removeWeakSpecies();
    
    // Create next generation
    std::vector<Genome> nextGen;
    nextGen.reserve(config.populationSize);
    
    // Elitism - copy best genomes
    auto eliteCount = static_cast<size_t>(config.populationSize * config.elitismRate);
    std::vector<size_t> indices(genomes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
        [&](size_t a, size_t b) { return genomes[a].getFitness() > genomes[b].getFitness(); });
    
    for (size_t i = 0; i < eliteCount && i < indices.size(); ++i) {
        nextGen.emplace_back(genomes[indices[i]]);
    }
    
    // Fill rest with offspring
    while (nextGen.size() < config.populationSize) {
        Genome child(selectParent().getConfig());
        
        // Crossover with second parent
        if (std::uniform_real_distribution<>(0, 1)(rng) < 0.75) {
            const auto& parent2 = selectParent();
            child = Genome::crossover(child, parent2, child.getConfig());
        }
        
        // Apply mutations
        child.mutateWeights();
        if (std::uniform_real_distribution<>(0, 1)(rng) < child.getConfig().newNodeRate) {
            child.addNodeMutation();
        }
        if (std::uniform_real_distribution<>(0, 1)(rng) < child.getConfig().newLinkRate) {
            child.addConnectionMutation();
        }
        
        nextGen.push_back(std::move(child));
    }
    
    genomes = std::move(nextGen);
}

const Genome& Population::getBestGenome() const {
    return *std::max_element(genomes.begin(), genomes.end(),
        [](const Genome& a, const Genome& b) { return a.getFitness() < b.getFitness(); });
}

const std::vector<Species>& Population::getSpecies() const {
    return species;
}

void Population::speciate() {
    species.clear();
    
    for (size_t i = 0; i < genomes.size(); ++i) {
        bool found = false;
        
        for (auto& spec : species) {
            if (Genome::compatibilityDistance(genomes[i], genomes[spec.getRepresentative()]) 
                < config.compatibilityThreshold) {
                spec.addGenome(i);
                found = true;
                break;
            }
        }
        
        if (!found) {
            Species newSpecies;
            newSpecies.setRepresentative(i);
            newSpecies.addGenome(i);
            species.push_back(std::move(newSpecies));
        }
    }
}

void Population::adjustFitness() {
    for (auto& spec : species) {
        size_t size = spec.getMembers().size();
        for (size_t idx : spec.getMembers()) {
            double adjustedFitness = genomes[idx].getFitness() / size;
            genomes[idx].setAdjustedFitness(adjustedFitness);
        }
    }
}

Genome& Population::selectParent() {
    std::vector<size_t> tournament;
    tournament.reserve(config.tournamentSize);
    
    std::uniform_int_distribution<size_t> dist(0, genomes.size() - 1);
    std::mt19937 rng(std::random_device{}());
    
    for (int32_t i = 0; i < config.tournamentSize; ++i) {
        tournament.push_back(dist(rng));
    }
    
    size_t winner = *std::max_element(tournament.begin(), tournament.end(),
        [this](size_t a, size_t b) { 
            return genomes[a].getAdjustedFitness() < genomes[b].getAdjustedFitness(); 
        });
        
    return genomes[winner];
}

int32_t Population::selectSpecies() {
    double totalFitness = 0.0;
    std::vector<double> speciesFitness;
    
    for (const auto& spec : species) {
        double fitness = 0.0;
        for (size_t idx : spec.getMembers()) {
            fitness += genomes[idx].getAdjustedFitness();
        }
        speciesFitness.push_back(fitness);
        totalFitness += fitness;
    }
    
    if (totalFitness <= 0.0) return -1;
    
    std::uniform_real_distribution<> dist(0, totalFitness);
    double spin = dist(rng);
    
    double sum = 0.0;
    for (size_t i = 0; i < species.size(); ++i) {
        sum += speciesFitness[i];
        if (sum > spin) return i;
    }
    
    return species.size() - 1;
}

void Population::removeStaleSpecies() {
    for (auto& spec : species) {
        double maxFitness = 0.0;
        for (size_t idx : spec.getMembers()) {
            maxFitness = std::max(maxFitness, genomes[idx].getFitness());
        }
        
        if (maxFitness > spec.getBestFitness()) {
            spec.setBestFitness(maxFitness);
            spec.resetStaleness();
        } else {
            spec.incrementStaleness();
        }
    }
    
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [](const Species& s) { return s.getStaleness() > 15; }),
        species.end());
}

void Population::removeWeakSpecies() {
    double totalAdjustedFitness = 0.0;
    for (const auto& spec : species) {
        for (size_t idx : spec.getMembers()) {
            totalAdjustedFitness += genomes[idx].getAdjustedFitness();
        }
    }
    
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [this, totalAdjustedFitness](const Species& s) {
                double speciesAdjustedFitness = 0.0;
                for (size_t idx : s.getMembers()) {
                    speciesAdjustedFitness += genomes[idx].getAdjustedFitness();
                }
                return (speciesAdjustedFitness / totalAdjustedFitness) * config.populationSize < 1.0;
            }),
        species.end());
}

}
}
