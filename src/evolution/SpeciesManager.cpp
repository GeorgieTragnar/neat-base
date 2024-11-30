#include "SpeciesManager.hpp"
#include "core/Population.hpp"
#include "core/Genome.hpp"
#include <algorithm>

namespace neat {
namespace evolution {

SpeciesManager::SpeciesManager(const Config& config)
    : config(config)
    , currentCompatibilityThreshold(config.compatibilityThreshold)
{}

void SpeciesManager::speciate(core::Population& population) {
    // Clear empty species
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [](const core::Species& s) { return s.getMembers().empty(); }),
        species.end()
    );
    
    // Assign each genome to a species
    for (size_t i = 0; i < population.size(); ++i) {
        assignGenomeToSpecies(population.getGenome(i), population);
    }
    
    adjustCompatibilityThreshold();
}

void SpeciesManager::updateSpecies(core::Population& population) {
    for (auto& s : species) {
        // Find max fitness in species
        double maxFitness = 0.0;
        for (auto idx : s.getMembers()) {
            maxFitness = std::max(maxFitness, population.getGenome(idx).getFitness());
        }
        s.updateFitness(maxFitness);
    }
}

void SpeciesManager::adjustFitness(core::Population& population) {
    for (auto& s : species) {
        double totalAdjustedFitness = 0.0;
        
        for (auto idx : s.getMembers()) {
            auto& genome = population.getGenome(idx);
            double adjustedFitness = genome.getFitness() / s.getMembers().size();
            genome.setAdjustedFitness(adjustedFitness);
            totalAdjustedFitness += adjustedFitness;
        }
        
        if (totalAdjustedFitness < config.minFitnessThreshold) {
            s.clearMembers();
        }
    }
}

void SpeciesManager::removeStagnantSpecies(core::Population& population) {
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [this](const core::Species& s) {
                return s.getStaleness() > config.stagnationThreshold &&
                       s.getBestFitness() < config.minFitnessThreshold;
            }),
        species.end()
    );
}

double SpeciesManager::getCompatibilityDistance(
    const core::Genome& genome1,
    const core::Genome& genome2) const {
    
    const auto& genes1 = genome1.getGenes();
    const auto& genes2 = genome2.getGenes();
    
    double weightDiff = 0.0;
    double actDiff = 0.0;
    int32_t matchingGenes = 0;
    int32_t disjointGenes = 0;
    
    size_t i = 0, j = 0;
    while (i < genes1.size() && j < genes2.size()) {
        if (genes1[i].innovation == genes2[j].innovation) {
            weightDiff += std::abs(genes1[i].weight - genes2[j].weight);

            // Activation difference
            actDiff += (genes1[i].activation.getType() != genes2[j].activation.getType()) ? 1.0 : 0.0;

            matchingGenes++;
            i++;
            j++;
        }
        else if (genes1[i].innovation < genes2[j].innovation) {
            disjointGenes++;
            i++;
        }
        else {
            disjointGenes++;
            j++;
        }
    }
    
    disjointGenes += (genes1.size() - i) + (genes2.size() - j);
    
    const double c1 = 1.0;  // Disjoint coefficient
    const double c2 = 1.0;  // Weight coefficient
    
    double distance = (c1 * disjointGenes) / std::max(genes1.size(), genes2.size());
    if (matchingGenes > 0) {
        distance += c2 * (weightDiff / matchingGenes);
        distance += config.activationDiffWeight * (actDiff / matchingGenes);
    }
    
    return distance;
}

void SpeciesManager::adjustCompatibilityThreshold() {
    if (species.size() < config.targetSpeciesCount) {
        currentCompatibilityThreshold *= 0.95;
    }
    else if (species.size() > config.targetSpeciesCount) {
        currentCompatibilityThreshold *= 1.05;
    }
}

void SpeciesManager::assignGenomeToSpecies(core::Genome& genome, core::Population& population) {
    static int32_t nextSpeciesId = 0;
    
    for (auto& s : species) {
        if (!s.getMembers().empty()) {
            const auto& representative = population.getGenome(s.getMembers()[0]);
            if (getCompatibilityDistance(genome, representative) < currentCompatibilityThreshold) {
                s.addGenome(genome.getSpecies());
                genome.setSpecies(nextSpeciesId);
                return;
            }
        }
    }
    
    // Create new species if no match found
    core::Species newSpecies;
    newSpecies.addGenome(genome.getSpecies());
    genome.setSpecies(nextSpeciesId);
    species.push_back(std::move(newSpecies));
    nextSpeciesId++;
}

const std::vector<core::Species>& SpeciesManager::getSpecies() const {
    return species;
}

}
}