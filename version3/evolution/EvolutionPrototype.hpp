#pragma once

#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <limits>
#include <cmath>
#include <type_traits>

// Core data structures
#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"

// Analysis and fitness
#include "version3/analysis/FitnessStrategy.hpp"
#include "version3/analysis/SpeciationControlUnit.hpp"

// Operators
#include "version3/operator/Init.hpp"
#include "version3/operator/Crossover.hpp"
#include "version3/operator/WeightMutation.hpp"
#include "version3/operator/ConnectionMutation.hpp"
#include "version3/operator/NodeMutation.hpp"
#include "version3/operator/ConnectionReactivation.hpp"
#include "version3/operator/CycleDetection.hpp"
#include "version3/operator/CompatibilityDistance.hpp"

// Population management
#include "version3/population/PopulationData.hpp"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/GenerationPlanner.hpp"
#include "version3/population/GenerationPlannerParams.hpp"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/population/ResolveParentIndices.hpp"
#include "version3/population/ReproductiveInstruction.hpp"

namespace Evolution {

// Simple results container
template<typename FitnessResultType>
class EvolutionResults {
public:
    EvolutionResults(
        Genome bestGenome,
        FitnessResultType bestFitness,
        std::vector<Genome> finalPopulation,
        std::vector<FitnessResultType> finalFitnessValues,
        uint32_t generationsCompleted
    ) : _bestGenome(std::move(bestGenome)),
        _bestFitness(bestFitness),
        _finalPopulation(std::move(finalPopulation)),
        _finalFitnessValues(std::move(finalFitnessValues)),
        _generationsCompleted(generationsCompleted) {}

    const Genome& getBestGenome() const { return _bestGenome; }
    FitnessResultType getBestFitness() const { return _bestFitness; }
    const std::vector<Genome>& getFinalPopulation() const { return _finalPopulation; }
    const std::vector<FitnessResultType>& getFinalFitnessValues() const { return _finalFitnessValues; }
    uint32_t getGenerationsCompleted() const { return _generationsCompleted; }

private:
    Genome _bestGenome;
    FitnessResultType _bestFitness;
    std::vector<Genome> _finalPopulation;
    std::vector<FitnessResultType> _finalFitnessValues;
    uint32_t _generationsCompleted;
};

// Main evolution prototype class
template<typename FitnessResultType>
class EvolutionPrototype {
public:
    EvolutionPrototype(
        std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
        uint32_t targetPopulationSize,
        const Operator::InitParams& initParams,
        const Population::GenerationPlannerParams& plannerParams,
        const Population::DynamicDataUpdateParams& updateParams,
        const Operator::CompatibilityDistanceParams& compatibilityParams,
        uint32_t randomSeed = std::random_device{}()
    );

    EvolutionResults<FitnessResultType> run(uint32_t maxGenerations);

private:
    // Core data containers
    std::array<std::vector<Genome>, 2> _population;
    std::array<std::multimap<FitnessResultType, Population::DynamicGenomeData>, 2> _genomeData;
    std::unordered_map<uint32_t, Population::DynamicSpeciesData> _speciesData;
    std::unordered_map<uint32_t, std::vector<Population::ReproductiveInstruction>> _instructionSets;
    std::shared_ptr<HistoryTracker> _historyTracker;

    // Strategy and parameters
    std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> _fitnessStrategy;
    uint32_t _targetPopulationSize;
    Operator::InitParams _initParams;
    Population::GenerationPlannerParams _plannerParams;
    Population::DynamicDataUpdateParams _updateParams;
    Operator::CompatibilityDistanceParams _compatibilityParams;
    std::mt19937 _rng;

    // No helper methods - everything in run()
};

// Template implementation
template<typename FitnessResultType>
EvolutionPrototype<FitnessResultType>::EvolutionPrototype(
    std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    uint32_t targetPopulationSize,
    const Operator::InitParams& initParams,
    const Population::GenerationPlannerParams& plannerParams,
    const Population::DynamicDataUpdateParams& updateParams,
    const Operator::CompatibilityDistanceParams& compatibilityParams,
    uint32_t randomSeed
) : _fitnessStrategy(std::move(fitnessStrategy)),
    _targetPopulationSize(targetPopulationSize),
    _initParams(initParams),
    _plannerParams(plannerParams),
    _updateParams(updateParams),
    _compatibilityParams(compatibilityParams),
    _rng(randomSeed),
    _historyTracker(std::make_shared<HistoryTracker>()) {
    
    // Basic parameter validation
    assert(_fitnessStrategy != nullptr);
    assert(_targetPopulationSize > 0);
    
    // Reserve capacity for initial population
    _population[0].reserve(_targetPopulationSize);
    
    // Create simple SpeciationControlUnit for fitness evaluation
    // TODO: Replace with proper implementation if fitness strategy needs speciation context
    class SimpleSpeciationControl : public Analysis::SpeciationControlUnit {
    public:
        std::vector<std::shared_ptr<const Analysis::Phenotype>> getChampions() const override { return {}; }
        std::shared_ptr<const Analysis::Phenotype> getBestChampion() const override { return nullptr; }
        std::shared_ptr<const Analysis::Phenotype> getRandomChampion() const override { return nullptr; }
        size_t getChampionCount() const override { return 0; }
    };
    SimpleSpeciationControl speciationControl;
    
    // Create initial genome
    Genome initialGenome = Operator::init(_historyTracker, _initParams);
    
    // Build phenotype needed for fitness evaluation
    Operator::phenotypeConstruct(initialGenome);
    
    const auto& cgenome = initialGenome;
    // Evaluate fitness
    FitnessResultType fitness = _fitnessStrategy->evaluate(
        cgenome.get_phenotype(), 
        speciationControl
    );
    
    // Determine species assignment
    uint32_t speciesId = Operator::compatibilityDistance(
        initialGenome, 
        _historyTracker, 
        _compatibilityParams
    );
    
    // Create genome metadata
    Population::DynamicGenomeData genomeData;
    genomeData.speciesId = speciesId;
    genomeData.protectionCounter = 0;
    genomeData.isUnderRepair = false;
    genomeData.isMarkedForElimination = false;
    
    // Add to population and genome data (generation 0)
    _population[0].push_back(std::move(initialGenome));
    _genomeData[0].insert({fitness, genomeData});
    // Important! the genome dynamic data must contain a valid pointer to its genome
    _genomeData[0].begin()->second.actualGenome = &(_population[0].back());
    
    // Run generation planner to create initial instruction sets
    _instructionSets = Population::generationPlanner(_speciesData, _plannerParams);

    // instructionSets should be empty because it runs on last generation data 
    // but that generation doesnt exist and the data is empty
    assert(_instructionSets.empty()); 

    Population::dynamicDataUpdate(_genomeData[0], _speciesData, _instructionSets, _updateParams);

    // at this point we should have bootstrapped the first evolution loop properly
}

template<typename FitnessResultType>
EvolutionResults<FitnessResultType> EvolutionPrototype<FitnessResultType>::run(uint32_t maxGenerations) {
    // Create a simple concrete implementation of SpeciationControlUnit
    class SimpleSpeciationControl : public Analysis::SpeciationControlUnit {
    public:
        std::vector<std::shared_ptr<const Analysis::Phenotype>> getChampions() const override {
            return {};
        }
        std::shared_ptr<const Analysis::Phenotype> getBestChampion() const override {
            return nullptr;
        }
        std::shared_ptr<const Analysis::Phenotype> getRandomChampion() const override {
            return nullptr;
        }
        size_t getChampionCount() const override {
            return 0;
        }
    };
    
    SimpleSpeciationControl speciationControl;

    for (uint32_t generation = 0; generation < maxGenerations; ++generation) {
        // // Convert _lastPopulation to vectors for indexed access
        // std::vector<Genome> parentGenomes;
        // std::vector<FitnessResultType> parentFitness;
        // for (const auto& [fitness, genome] : _lastPopulation) {
        //     parentGenomes.push_back(genome);
        //     parentFitness.push_back(fitness);
        // }

        // // Clear current population for new generation
        // _currentPopulation.clear();

        // // Iterate through instruction sets and execute
        // for (const auto& [speciesId, instructions] : _instructionSets) {
        //     for (const auto& instruction : instructions) {
        //         // if instruction set for speciesId does not exist that species is marked for death
        //         if (instruction.globalParentIndices.empty()) continue;
                
        //         // Get first parent
        //         size_t parentIdx1 = instruction.globalParentIndices[0];
        //         if (parentIdx1 >= parentGenomes.size()) continue;
                
        //         Genome offspring = parentGenomes[parentIdx1]; // Copy first parent
                
        //         switch (instruction.operationType) {
        //             case Population::OperationType::PRESERVE:
        //                 // Do nothing - offspring is already a copy of parent
        //                 break;
                        
        //             case Population::OperationType::MUTATE_UNPROTECTED:
        //             case Population::OperationType::MUTATE_PROTECTED: {
        //                 // Apply weight mutation using instruction parameters
        //                 Operator::WeightMutationParams weightParams(
        //                     instruction.evolutionParams.mutationRate,
        //                     0.1, // replacement rate
        //                     0.5, // perturbation strength  
        //                     2.0, // weight range
        //                     Operator::WeightMutationParams::MutationType::MIXED
        //                 );
        //                 offspring = Operator::weightMutation(offspring, weightParams);
        //                 break;
        //             }
                    
        //             case Population::OperationType::CROSSOVER: {
        //                 if (instruction.globalParentIndices.size() >= 2) {
        //                     size_t parentIdx2 = instruction.globalParentIndices[1];
        //                     if (parentIdx2 < parentGenomes.size()) {
        //                         Operator::CrossoverParams crossoverParams;
        //                         offspring = Operator::crossover(
        //                             offspring, parentFitness[parentIdx1],
        //                             parentGenomes[parentIdx2], parentFitness[parentIdx2],
        //                             crossoverParams
        //                         );
        //                     }
        //                 }
        //                 break;
        //             }
        //         }
                
        //         // Add offspring to current population (fitness will be evaluated later)
        //         _currentPopulation.insert({FitnessResultType{}, offspring});
        //     }
        //     // end of instruction set loop
        // }

        // // end of generation loop
    }
    
    // // Build results directly
    // if (_currentPopulation.empty()) {
    //     // Create a dummy genome for empty population case
    //     Genome dummyGenome = Operator::init(_historyTracker, _initParams);
    //     return EvolutionResults<FitnessResultType>(
    //         std::move(dummyGenome),
    //         FitnessResultType{},
    //         std::vector<Genome>{},
    //         std::vector<FitnessResultType>{},
    //         maxGenerations
    //     );
    // }
    
    // // Extract best genome (last in multimap - highest fitness in reverse order)
    // auto bestIt = _currentPopulation.rbegin();
    // Genome bestGenome = bestIt->second;
    // FitnessResultType bestFitness = bestIt->first;
    
    // // Extract final population
    // std::vector<Genome> finalPopulation;
    // std::vector<FitnessResultType> finalFitnessValues;
    
    // for (const auto& [fitness, genome] : _currentPopulation) {
    //     finalPopulation.push_back(genome);
    //     finalFitnessValues.push_back(fitness);
    // }
    
    // return EvolutionResults<FitnessResultType>(
    //     std::move(bestGenome),
    //     bestFitness,
    //     std::move(finalPopulation),
    //     std::move(finalFitnessValues),
    //     maxGenerations
    // );

    Genome dummyGenome = Operator::init(_historyTracker, _initParams);
    return EvolutionResults<FitnessResultType>(
        std::move(dummyGenome),
        FitnessResultType{},
        std::vector<Genome>{},
        std::vector<FitnessResultType>{},
        maxGenerations
    );
}



} // namespace Evolution