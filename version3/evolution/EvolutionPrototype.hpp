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
    std::array<std::unordered_map<uint32_t, std::vector<Population::ReproductiveInstruction>>, 2> _instructionSets;
    std::shared_ptr<HistoryTracker> _historyTracker;

    uint32_t _lastGeneration;
    uint32_t _currentGeneration;

    // Strategy and parameters
    std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> _fitnessStrategy;
    uint32_t _targetPopulationSize;
    Operator::InitParams _initParams;
    Population::GenerationPlannerParams _plannerParams;
    Population::DynamicDataUpdateParams _updateParams;
    Operator::CompatibilityDistanceParams _compatibilityParams;
    std::mt19937 _rng;

protected:
    Genome createOffspring(
        const Population::ReproductiveInstruction& instruction,
        const std::vector<decltype(_genomeData[_lastGeneration].begin())>& indexToIterator,
        bool& hasCycles);
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
    genomeData.genomeIndex = 0;
    
    // Add to population and genome data (generation 0)
    _population[0].push_back(std::move(initialGenome));
    _genomeData[0].insert({fitness, genomeData});
    
    // Run generation planner to create initial instruction sets
    _instructionSets[0] = Population::generationPlanner(_speciesData, _plannerParams);

    // instructionSets should be empty because it runs on last generation data 
    // but that generation doesnt exist and the data is empty
    assert(_instructionSets[0].empty()); 

    Population::dynamicDataUpdate(_genomeData[0], _speciesData, _instructionSets[0], _updateParams);

    _currentGeneration = 0;
    // at this point we should have bootstrapped the first evolution loop properly
}

template<typename FitnessResultType>
Genome EvolutionPrototype<FitnessResultType>::createOffspring(
    const Population::ReproductiveInstruction& instruction,
    const std::vector<decltype(_genomeData[_lastGeneration].begin())>& indexToIterator,
    bool& hasCycles) 
{
    switch (instruction.operationType) {
        case Population::OperationType::PRESERVE: {
            return _population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex];
        }
        
        case Population::OperationType::MUTATE_UNPROTECTED:
        case Population::OperationType::MUTATE_PROTECTED: {                        
            // Randomly select mutation type
            // TODO: currently not allowing connection reactivation
            // we are missing how to check if there deactivated connections
            std::uniform_int_distribution<int> mutationDist(0, 2);
            int mutationType = mutationDist(_rng);
            
            switch (mutationType) {
                case 0: { // Weight mutation
                    Operator::WeightMutationParams weightParams(
                        instruction.evolutionParams.mutationRate,
                        0.1, 0.5, 2.0, Operator::WeightMutationParams::MutationType::MIXED
                    );
                    Genome offspring = Operator::weightMutation(_population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex]
                        , weightParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        Operator::phenotypeUpdateWeight(offspring);
                    }
                    return std::move(offspring);
                }
                
                case 1: { // Node mutation
                    NodeGeneAttributes nodeAttribs{ActivationType::SIGMOID}; // TODO: Configure
                    Operator::NodeMutationParams nodeParams(nodeAttribs);
                    Genome offspring = Operator::nodeMutation(_population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex]
                        , _historyTracker, nodeParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        Operator::phenotypeUpdateNode(offspring);
                    }
                    return std::move(offspring);
                }
                
                case 2: { // Connection mutation
                    Operator::ConnectionMutationParams connParams(
                        instruction.evolutionParams.mutationRate, 2.0,
                        Operator::ConnectionMutationParams::NetworkTopology::FEED_FORWARD
                    );
                    Genome offspring = Operator::connectionMutation(_population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex]
                        , _historyTracker, connParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        Operator::phenotypeUpdateConnection(offspring);
                    }
                    return std::move(offspring);
                }
                
                case 3: { // Connection reactivation
                    // TODO: currently disabled
                    // TODO: Check if disabled connections exist first
                    Operator::ConnectionReactivationParams reactParams(
                        Operator::ConnectionReactivationParams::SelectionStrategy::RANDOM
                    );
                    Genome offspring = Operator::connectionReactivation(_population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex]
                        , reactParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        Operator::phenotypeUpdateConnection(offspring);
                    }
                    return std::move(offspring);
                }
                default: {
                    assert(false && "we cannot do nothing for an instruction set");
                    return Operator::init(_historyTracker, _initParams);;
                }
            }
        }
        
        case Population::OperationType::CROSSOVER: {                        
            Operator::CrossoverParams crossoverParams;
            Genome offspring = Operator::crossover(
                _population[_lastGeneration][indexToIterator[instruction.globalParentIndices[0]]->second.genomeIndex],
                indexToIterator[instruction.globalParentIndices[0]]->first,
                _population[_lastGeneration][indexToIterator[instruction.globalParentIndices[1]]->second.genomeIndex],
                indexToIterator[instruction.globalParentIndices[1]]->first,
                crossoverParams
            );

            Operator::CycleDetectionParams cycleParams;
            hasCycles = Operator::hasCycles(offspring, cycleParams);

            if (!hasCycles) {
                Operator::phenotypeConstruct(offspring);
            }
            return std::move(offspring);
        }

        default: {
            assert(false && "we cannot do nothing for an instruction set");
            return Operator::init(_historyTracker, _initParams);;
        }
    }
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

        _lastGeneration = _currentGeneration;
        _currentGeneration = (_currentGeneration + 1) % 2;

        // Firstly calculate how many genomes the current generation is going to have
        // Which is all instruction sets plus all species that dont have instructions sets
        // but have population size in dynamic species data

        // then we reserve this size inside the genome vector we are going to use
        // to prevent reordering for performance, the ordering of genomes 
        // is important because these indices are stored inside the dynamic genome data

        // afterwards we start executing the instruction sets and fill the population up
        // however we need to check for cycles in order to get fitness result
        // and create new entries in the multimap of fitness results for dynamic genome data
        // we also need to get compatibility distance result before we fully fillup 
        // the dynamic genome data for inserting into the multimap
        // we need to keep in mind that we have to work on the copy of parent dynamic genome data
        // except when the instruction set was crossover, because that indicates we should
        // create empty dynamic genome data

        // afterwards we do all of the population operations
        // which is generationplanner, speciesgrouping, resolveparentindices and dynamicdataupdate

        // this makes the loop complete as the dynamic data update is the loops finishing operation


        // Calculate total genomes for current generation
        uint32_t totalGenomesThisGeneration = 0;

        for (const auto& [speciesId, speciesData] : _speciesData) {
            auto instructionIt = _instructionSets[_lastGeneration].find(speciesId);
            if (instructionIt != _instructionSets[_lastGeneration].end()) {
                // Species has instruction sets (active or marked for elimination)
                totalGenomesThisGeneration += static_cast<uint32_t>(instructionIt->second.size());
            } else {
                // Species not in instruction sets = new/rediscovered, carry forward
                totalGenomesThisGeneration += speciesData.currentPopulationSize;
            }
        }
        _population[_currentGeneration].clear();
        _population[_currentGeneration].reserve(totalGenomesThisGeneration);

        // Build index-to-iterator mapping once per generation
        std::vector<decltype(_genomeData[_lastGeneration].begin())> indexToIterator;
        indexToIterator.reserve(_genomeData[_lastGeneration].size());
        for (auto it = _genomeData[_lastGeneration].begin(); it != _genomeData[_lastGeneration].end(); ++it) {
            indexToIterator.push_back(it);
        }
        
        // Execute instruction sets and create new generation
        for (const auto& [speciesId, instructions] : _instructionSets[_lastGeneration]) {
            for (const auto& instruction : instructions) {
                // Execute the instruction to create offspring
                bool hasCycles = false;
                Genome offspring = createOffspring(instruction, indexToIterator, hasCycles);

                // Create dynamic genome data for offspring
                Population::DynamicGenomeData genomeData = indexToIterator[instruction.globalParentIndices[0]]->second;
                genomeData.genomeIndex = static_cast<uint32_t>(_population[_currentGeneration].size());
                
                if (instruction.operationType == Population::OperationType::CROSSOVER) {
                    // Fresh data for crossover offspring
                    genomeData.protectionCounter = 0;
                    genomeData.isMarkedForElimination = false;
                }
                
                if (!hasCycles) {
                    // Evaluate fitness
                    const auto& coffspring = offspring;
                    FitnessResultType fitness = _fitnessStrategy->evaluate(
                        coffspring.get_phenotype(), speciationControl
                    );
                    
                    // Get species assignment
                    uint32_t offspringSpeciesId = Operator::compatibilityDistance(
                        offspring, _historyTracker, _compatibilityParams
                    );

                    genomeData.speciesId = offspringSpeciesId;
                    genomeData.isUnderRepair = false;

                    _genomeData[_currentGeneration].insert({fitness, genomeData});
                }
                else {
                    // if the offspring needs repair and cant be evaluated
                    // then assign it the species of the first parent
                    // even for crossover
                    genomeData.isUnderRepair = true;
                    _genomeData[_currentGeneration].insert({FitnessResultType(), genomeData});
                }
                
                // Add to current generation
                _population[_currentGeneration].push_back(std::move(offspring));
            }
            // instruction set end of loop
        }

        // Population operations
        // 1. Generate instruction sets for next generation based on current species data
        _instructionSets[_currentGeneration] = Population::generationPlanner(_speciesData, _plannerParams);

        // 2. Group current generation by species to get index mappings
        auto speciesGrouping = Population::speciesGrouping(_genomeData[_currentGeneration], _speciesData);

        // 3. Resolve relative parent indices to global indices in instruction sets
        for (auto& [speciesId, instructions] : _instructionSets[_currentGeneration]) {
            auto groupingIt = speciesGrouping.find(speciesId);
            // if species grouping doesnt contain these speciesId then there is a bug
            assert(groupingIt != speciesGrouping.end());
            Population::resolveParentIndices(instructions, groupingIt->second);
        }

        // 4. Update dynamic data for current generation (finishing operation of evolution loop)
        Population::dynamicDataUpdate(_genomeData[_currentGeneration], _speciesData, _instructionSets[_currentGeneration], _updateParams);
        // end of generation loop
    }
    
    // Build results directly
    if (_genomeData[_currentGeneration].empty()) {
        // Create a dummy genome for empty population case
        Genome dummyGenome = Operator::init(_historyTracker, _initParams);
        return EvolutionResults<FitnessResultType>(
            std::move(dummyGenome),
            FitnessResultType{},
            std::vector<Genome>{},
            std::vector<FitnessResultType>{},
            maxGenerations
        );
    }
    
    // Extract best genome (last in multimap - highest fitness in reverse order)
    auto bestIt = _genomeData[_currentGeneration].rbegin();
    Genome bestGenome = _population[_currentGeneration][bestIt->second.genomeIndex];
    FitnessResultType bestFitness = bestIt->first;
    
    // Extract final population
    std::vector<Genome> finalPopulation;
    std::vector<FitnessResultType> finalFitnessValues;

    std::vector<decltype(_genomeData[_currentGeneration].begin())> indexToIterator;
    indexToIterator.reserve(_genomeData[_currentGeneration].size());
    for (auto it = _genomeData[_currentGeneration].begin(); it != _genomeData[_currentGeneration].end(); ++it) {
        indexToIterator.push_back(it);
    }
    
    for (const auto& [fitness, genomeData] : _genomeData[_currentGeneration]) {
        finalPopulation.push_back(_population[_currentGeneration][genomeData.genomeIndex]);
        finalFitnessValues.push_back(fitness);
    }
    
    return EvolutionResults<FitnessResultType>(
        std::move(bestGenome),
        bestFitness,
        std::move(finalPopulation),
        std::move(finalFitnessValues),
        maxGenerations
    );
}



} // namespace Evolution