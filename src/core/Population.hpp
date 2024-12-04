// population.hpp
#pragma once
#include <cstdint>
#include <functional>
#include <random>
#include "core/Genome.hpp"
#include "evolution/MutationOperator.hpp"
#include "evolution/CrossoverOperator.hpp"

namespace neat {
namespace core {

class Genome;
class Species;

class Population {
public:
    struct Config {
        // Core population parameters
        int32_t inputSize;
        int32_t outputSize;
        int32_t populationSize;
        double survivalThreshold;
        double compatibilityThreshold;
        int32_t speciesTargetSize;
        int32_t tournamentSize;
        double elitismRate;

        ActivationGene::Config activationConfig;

        // Factory methods for child configs
        Genome::Config createGenomeConfig() const {
            return Genome::Config{
                inputSize,
                outputSize,
                0.8,    // Weight mutation rate
                0.9,    // Weight perturbation rate
                0.03,   // New node rate
                0.05,   // New link rate
                2.0,    // Weight perturbation range
                5.0,    // New weight range
                5,      // Max hidden nodes
                activationConfig // Pass through activation config
            };
        }

        evolution::MutationOperator::Config createMutationConfig() const {
            return evolution::MutationOperator::Config(
                0.8,    // Weight mutation rate
                0.9,    // Weight perturbation rate
                0.03,   // New node rate
                0.05,   // New connection rate
                2.0,    // Perturbation range
                5.0,    // New weight range
                activationConfig.mutationRate,
                activationConfig.interTierMutationRate,
                activationConfig,  // Pass through activation config
                false   // Allow recurrent
            );
        }

        evolution::CrossoverOperator::Config createCrossoverConfig() const {
            return evolution::CrossoverOperator::Config(
                0.7,    // matchingGeneInheritanceRate - favor inherited genes
                true,   // inheritDisabledGenes - allow inheritance of disabled genes
                0.25,   // disabledGeneReenableRate - conservative re-enabling
                0.6     // matchingActivationInheritanceRate - favor matched activations
            );
        }

    private:
        // Helper to create activation config used by mutation
        core::ActivationGene::Config createActivationConfig() const {
            return core::ActivationGene::Config{
                0.1,    // Mutation rate
                0.2,    // Inter-tier rate
                0.8,    // Basic tier prob
                0.15,   // Advanced tier prob
                1.0     // Compatibility weight
            };
        }

    public:
        Config() = delete;
        Config(int32_t inputs, int32_t outputs, 
               int32_t popSize = 500,
               double survivalThresh = 0.2,
               double compatThresh = 3.0,
               int32_t speciesTarget = 10,
               int32_t tournSize = 3,
               double elitismR = 0.1,
               const ActivationGene::Config& actConfig = ActivationGene::Config(
                   0.1,    // mutation rate
                   0.2,    // inter tier rate
                   0.8,    // basic tier prob
                   0.15,   // advanced tier prob
                   1.0     // compatibility weight
               ))
            : inputSize(inputs)
            , outputSize(outputs)
            , populationSize(popSize)
            , survivalThreshold(survivalThresh)
            , compatibilityThreshold(compatThresh)
            , speciesTargetSize(speciesTarget)
            , tournamentSize(tournSize)
            , elitismRate(elitismR) 
            , activationConfig(actConfig) {}
    };
    
    Population(const Config& config);
    
    void evolve(const std::function<double(Genome&)>& fitnessFunction);
    const Genome& getBestGenome() const;
    const std::vector<Species>& getSpecies() const;
    
    const Genome& getGenome(size_t idx) const { return genomes[idx]; }
    Genome& getGenome(size_t idx) { return genomes[idx]; }
    
    size_t size() const { return genomes.size(); }
    const Config& getConfig() const { return config; }
    
private:
    std::vector<Genome> genomes;
    std::vector<Species> species;
    Config config;
    
    void speciate();
    void adjustFitness();
    Genome& selectParent();
    int32_t selectSpecies();
    void removeStaleSpecies();
    void removeWeakSpecies();

    evolution::MutationOperator mutationOp;
    std::mt19937 rng;
};

}
}
