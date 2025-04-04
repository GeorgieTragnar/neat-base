#include "Population.hpp"
#include "Genome.hpp"
#include "Species.hpp"
#include "evolution/CrossoverOperator.hpp"
#include <algorithm>
#include <random>
#include <iostream>

#include "logger/Logger.hpp"
static auto logger = LOGGER("core::Population");

namespace neat {
namespace core {

Population::Population(int32_t inputSize, int32_t outputSize, const Config& config)
    : config(config)
    , nextSpeciesId(0)
    , rng(std::random_device{}()) {
    
    LOG_DEBUG("Creating population with {} inputs and {} outputs", inputSize, outputSize);
              
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
        
        LOG_DEBUG("Genome {} initialized with {} inputs and {} outputs", i, inputCount, outputCount);
                  
        genomes.push_back(std::move(genome));
    }
    
    LOG_DEBUG("Population initialized with {} genomes", genomes.size());
              
    speciate();
}

void Population::evolve(const std::function<double(const Genome&)>& fitnessFunction) {
    LOG_TRACE("Starting evolution cycle");

    // Create operators
    evolution::CrossoverOperator crossover(evolution::CrossoverOperator::Config{});

    // Evaluate fitness
    for (size_t i = 0; i < genomes.size(); i++) {
        auto& genome = genomes[i];
        genome.validate();  // Validate before fitness evaluation
        double fitness = fitnessFunction(genome);
        genome.setFitness(fitness);
    }
    
    // Update species and adjust fitness
    LOG_TRACE("Updating species...");
    speciate();
    adjustFitness();
    removeStaleSpecies();
    removeWeakSpecies();
    
    // Create next generation
    std::vector<Genome> nextGen;
    nextGen.reserve(config.populationSize);
    
    // Elitism - copy best genomes
    auto eliteCount = static_cast<size_t>(config.populationSize * config.elitismRate);
    LOG_TRACE("Selecting {} elites", eliteCount);
    
    std::vector<size_t> indices(genomes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
        [&](size_t a, size_t b) { return genomes[a].getFitness() > genomes[b].getFitness(); });
    
    for (size_t i = 0; i < eliteCount && i < indices.size(); ++i) {
        nextGen.emplace_back(genomes[indices[i]]);
        nextGen.back().validate();  // Verify elite copy
    }
    
    // Fill rest with offspring
    while (nextGen.size() < config.populationSize) {
        // Select and copy first parent
        const auto& parent1 = selectParent();
        Genome child(parent1);  // Use copy constructor to get complete genome
        child.validate();  // Verify copy
        
        // Crossover with second parent
        if (std::uniform_real_distribution<>(0, 1)(rng) < 0.75) {
            const auto& parent2 = selectParent();
            child = crossover.crossover(parent1, parent2);
            child.validate();  // Verify after crossover
        }
        
        LOG_TRACE("Before mutations - Child has {} inputs and {} outputs",
            std::count_if(child.getNodes().begin(), child.getNodes().end(),
                      [](const auto& pair) { return pair.second == core::ENodeType::INPUT; }),
            std::count_if(child.getNodes().begin(), child.getNodes().end(),
                      [](const auto& pair) { return pair.second == core::ENodeType::OUTPUT; }));
        
        // Apply mutations
        child.mutateWeights();
        if (std::uniform_real_distribution<>(0, 1)(rng) < child.getConfig().newNodeRate) {
            LOG_TRACE("Attempting node mutation...");
            child.addNodeMutation();
        }
        if (std::uniform_real_distribution<>(0, 1)(rng) < child.getConfig().newLinkRate) {
            LOG_TRACE("Attempting connection mutation...");
            child.addConnectionMutation();
        }
        
        LOG_TRACE("After mutations - Child has {} inputs and {} outputs",
            std::count_if(child.getNodes().begin(), child.getNodes().end(),
                      [](const auto& pair) { return pair.second == core::ENodeType::INPUT; }),
            std::count_if(child.getNodes().begin(), child.getNodes().end(),
                      [](const auto& pair) { return pair.second == core::ENodeType::OUTPUT; }));
        
        child.validate();  // Final validation before adding to next generation
        nextGen.push_back(std::move(child));
    }
    
    genomes = std::move(nextGen);
    
    // Validate final population
    for (auto& genome : genomes) {
        genome.validate();
    }

    selectSpeciesRepresentatives();
    
    LOG_TRACE("Evolution cycle complete");
}

const Genome& Population::getBestGenome() const {
    return *std::max_element(genomes.begin(), genomes.end(),
        [](const Genome& a, const Genome& b) { return a.getFitness() < b.getFitness(); });
}

const std::vector<Species>& Population::getSpecies() const {
    return species;
}

void Population::speciate() {
    // Clear all species members but keep the species
    for (auto& spec : species) {
        spec.clearMembers();
        // Don't reset representatives yet - we need them for comparison
    }
    
    // Mark all genomes as unassigned initially
    for (size_t i = 0; i < genomes.size(); ++i) {
        genomes[i].setSpecies(-1);
    }
    
    // First pass: assign genomes to existing species using representatives
    for (size_t i = 0; i < genomes.size(); ++i) {
        bool found = false;
        
        for (size_t s = 0; s < species.size(); ++s) {
            // Only consider species that have a representative
            if (species[s].hasRepresentative()) {
                if (Genome::compatibilityDistance(
                        genomes[i], 
                        genomes[species[s].getRepresentative()], 
                        config.genomeConfig) < config.compatibilityThreshold) {
                    
                    species[s].addGenome(i);
                    genomes[i].setSpecies(species[s].getId());
                    found = true;
                    break;
                }
            }
        }
        
        // If not found in any existing species, create a new one
        if (!found) {
            Species newSpecies;
            int32_t newId = nextSpeciesId++;
            
            newSpecies.addGenome(i);
            newSpecies.setRepresentative(i);  // Set this genome as the representative
            newSpecies.setId(newId);
            genomes[i].setSpecies(newId);
            
            species.push_back(std::move(newSpecies));
        }
    }
    
    // Remove empty species
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [](const Species& s) { return s.getMembers().empty(); }),
        species.end());
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

void Population::selectSpeciesRepresentatives() {
    for (auto& spec : species) {
        if (!spec.getMembers().empty()) {
            // Select a random member as representative for the next generation
            std::uniform_int_distribution<size_t> dist(0, spec.getMembers().size() - 1);
            size_t randomIdx = dist(rng);
            size_t representativeGenomeIdx = spec.getMembers()[randomIdx];
            spec.setRepresentative(representativeGenomeIdx);
        } else {
            // No members, so no representative
            spec.resetRepresentative();
        }
    }
}

}
}
