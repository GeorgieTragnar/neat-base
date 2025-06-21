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
#include "version3/data/NodeGene.hpp"
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
#include "version3/operator/RepairOperator.hpp"
#include "version3/operator/EmptyDeltas.hpp"
#include "version3/operator/HasDisabledConnections.hpp"

// Population management
#include "version3/population/PopulationData.hpp"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/population/PlotElites.hpp"
#include "version3/population/PlotCrossover.hpp"

#include "../logger/Logger.hpp"

namespace Evolution {

auto logger = LOGGER("evolution.EvolutionPrototype");

// Mutation probability parameters
struct MutationProbabilityParams {
    double weightMutationProbability;        // High probability for weight mutations
    double nodeMutationProbability;         // Low probability for structural changes
    double connectionMutationProbability;   // Low probability for structural changes
    double connectionReactivationProbability; // Very low probability for reactivation
    
    MutationProbabilityParams() = delete;
    
    MutationProbabilityParams(double weightProb, double nodeProb, double connProb, double reactProb)
        : weightMutationProbability(weightProb), 
          nodeMutationProbability(nodeProb),
          connectionMutationProbability(connProb),
          connectionReactivationProbability(reactProb) {
        
        // Normalize probabilities to sum to 1.0
        double total = weightProb + nodeProb + connProb + reactProb;
        if (total > 0) {
            weightMutationProbability /= total;
            nodeMutationProbability /= total;
            connectionMutationProbability /= total;
            connectionReactivationProbability /= total;
        }
    }
};
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
        const Population::DynamicDataUpdateParams& updateParams,
        const Population::PlotElitesParams& eliteParams,
        const Population::PlotCrossoverParams& crossoverParams,
        const Operator::CompatibilityDistanceParams& compatibilityParams,
        const Operator::RepairOperatorParams& repairParams,
        const MutationProbabilityParams& mutationParams,
        const Operator::WeightMutationParams& weightMutationParams,
        const Operator::CrossoverParams& crossoverOperatorParams,
        const Operator::CycleDetectionParams& cycleDetectionParams,
        const Operator::NodeMutationParams& nodeMutationParams,
        const Operator::ConnectionMutationParams& connectionMutationParams,
        const Operator::ConnectionReactivationParams& connectionReactivationParams,
        uint32_t randomSeed = std::random_device{}()
    );

    EvolutionResults<FitnessResultType> run(uint32_t maxGenerations);

private:

    // Core data containers - simplified architecture
    std::array<std::vector<Genome>, 2> _population;
    std::array<std::vector<Population::DynamicGenomeData>, 2> _genomeData;
    std::array<std::multimap<FitnessResultType, size_t>, 2> _fitnessResults;
    std::unordered_map<uint32_t, Population::DynamicSpeciesData> _speciesData;
    std::shared_ptr<HistoryTracker> _historyTracker;

    uint32_t _lastGeneration;
    uint32_t _currentGeneration;

    // Strategy and parameters
    std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> _fitnessStrategy;
    uint32_t _targetPopulationSize;
    Operator::InitParams _initParams;
    Population::DynamicDataUpdateParams _updateParams;
    Population::PlotElitesParams _eliteParams;
    Population::PlotCrossoverParams _crossoverParams;
    Operator::CompatibilityDistanceParams _compatibilityParams;
    Operator::RepairOperatorParams _repairParams;
    MutationProbabilityParams _mutationParams;
    Operator::WeightMutationParams _weightMutationParams;
    Operator::CrossoverParams _crossoverOperatorParams;
    Operator::CycleDetectionParams _cycleDetectionParams;
    Operator::NodeMutationParams _nodeMutationParams;
    Operator::ConnectionMutationParams _connectionMutationParams;
    Operator::ConnectionReactivationParams _connectionReactivationParams;
    std::mt19937 _rng;

};

// Template implementation
template<typename FitnessResultType>
EvolutionPrototype<FitnessResultType>::EvolutionPrototype(
    std::unique_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    uint32_t targetPopulationSize,
    const Operator::InitParams& initParams,
    const Population::DynamicDataUpdateParams& updateParams,
    const Population::PlotElitesParams& eliteParams,
    const Population::PlotCrossoverParams& crossoverParams,
    const Operator::CompatibilityDistanceParams& compatibilityParams,
    const Operator::RepairOperatorParams& repairParams,
    const MutationProbabilityParams& mutationParams,
    const Operator::WeightMutationParams& weightMutationParams,
    const Operator::CrossoverParams& crossoverOperatorParams,
    const Operator::CycleDetectionParams& cycleDetectionParams,
    const Operator::NodeMutationParams& nodeMutationParams,
    const Operator::ConnectionMutationParams& connectionMutationParams,
    const Operator::ConnectionReactivationParams& connectionReactivationParams,
    uint32_t randomSeed
) : _fitnessStrategy(std::move(fitnessStrategy)),
    _targetPopulationSize(targetPopulationSize),
    _initParams(initParams),
    _updateParams(updateParams),
    _eliteParams(eliteParams),
    _crossoverParams(crossoverParams),
    _compatibilityParams(compatibilityParams),
    _repairParams(repairParams),
    _mutationParams(mutationParams),
    _weightMutationParams(weightMutationParams),
    _crossoverOperatorParams(crossoverOperatorParams),
    _cycleDetectionParams(cycleDetectionParams),
    _nodeMutationParams(nodeMutationParams),
    _connectionMutationParams(connectionMutationParams),
    _connectionReactivationParams(connectionReactivationParams),
    _rng(randomSeed),
    _historyTracker(std::make_shared<HistoryTracker>()) {
    
    // Basic parameter validation
    assert(_fitnessStrategy != nullptr);
    assert(_targetPopulationSize > 0);
    
    // Reserve capacity for initial population
    _population[0].reserve(_targetPopulationSize);
    _population[1].reserve(_targetPopulationSize);
    _genomeData[0].reserve(_targetPopulationSize);
    _genomeData[1].reserve(_targetPopulationSize);
    
    // Create initial genome
    Genome initialGenome = Operator::init(_historyTracker, _initParams);
    Operator::phenotypeConstruct(initialGenome);
    
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
    _genomeData[0].push_back(genomeData);
    
    _currentGeneration = 0;
    _lastGeneration = 1;
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

        // Clear current generation containers
        _population[_currentGeneration].clear();
        _genomeData[_currentGeneration].clear();
        _fitnessResults[_currentGeneration].clear();

        const size_t populationSize = _population[_lastGeneration].size();
        
        // Reserve capacity for growth
        if (_population[_currentGeneration].capacity() < populationSize * 2) {
            _population[_currentGeneration].reserve(populationSize * 2);
            _genomeData[_currentGeneration].reserve(populationSize * 2);
        }

        LOG_DEBUG("Generation {}: Starting with {} genomes", generation, populationSize);

        // Phase 1: 1:1 Evolution - Copy and mutate each genome from last generation
        for (size_t i = 0; i < populationSize; ++i) {
            const Genome& parentGenome = _population[_lastGeneration][i];
            const Population::DynamicGenomeData& parentData = _genomeData[_lastGeneration][i];
            
            if (parentData.isMarkedForElimination) {
                _genomeData[_currentGeneration][i] = parentData;
                continue;
            }
            
            // Create corresponding genome data
            Population::DynamicGenomeData offspringData = parentData;
            
            bool hasCycles = false;
            
            // Copy parent genome
            Genome offspring = parentGenome;
            
            // Check if parent is under repair - if so, attempt repair instead of evolution
            if (parentData.isUnderRepair) {
                offspring = Operator::repair(offspring, offspringData, _repairParams);
                
                // Check if repair was successful
                hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                
                if (!hasCycles) {
                    // Repair successful - no longer under repair
                    offspringData.isUnderRepair = false;
                } else {
                    // Repair failed - still under repair
                    offspringData.isUnderRepair = true;
                }
            } else {
                // Normal evolution path - apply mutations based on probability
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                double random = dist(_rng);
                
                if (random < _mutationParams.weightMutationProbability) {
                    // Weight mutation
                    offspring = Operator::weightMutation(offspring, _weightMutationParams);
                    hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                    if (!hasCycles) {
                        Operator::phenotypeUpdateWeight(offspring);
                    }
                } else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability) {
                    // Node mutation
                    offspring = Operator::nodeMutation(offspring, _historyTracker, _nodeMutationParams);
                    hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                    if (!hasCycles) {
                        Operator::phenotypeUpdateNode(offspring);
                    }
                } else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability + _mutationParams.connectionMutationProbability) {
                    // Connection mutation
                    offspring = Operator::connectionMutation(offspring, _historyTracker, _connectionMutationParams);
                    hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                    if (!hasCycles) {
                        Operator::phenotypeUpdateConnection(offspring);
                    }
                } else {
                    // Connection reactivation
                    // Check if genome has any disabled connections using analytical operator
                    
                    if (Operator::hasDisabledConnections(offspring)) {
                        offspring = Operator::connectionReactivation(offspring, _connectionReactivationParams);
                        hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        if (!hasCycles) {
                            Operator::phenotypeUpdateConnection(offspring);
                        }
                    } else {
                        // Fallback to weight mutation
                        offspring = Operator::weightMutation(offspring, _weightMutationParams);
                        hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        if (!hasCycles) {
                            Operator::phenotypeUpdateWeight(offspring);
                        }
                    }
                }
            }
            
            // Mark as under repair if cycles detected
            if (hasCycles) {
                offspringData.isUnderRepair = true;
            } else {
                offspringData.isUnderRepair = false;
            }
            
            // Add to current generation
            _population[_currentGeneration].push_back(std::move(offspring));
            _genomeData[_currentGeneration].push_back(offspringData);
        }

        // Phase 2: Fitness Evaluation - Build fitness results multimap
        for (size_t i = 0; i < _population[_currentGeneration].size(); ++i) {
            const Genome& genome = _population[_currentGeneration][i];
            const Population::DynamicGenomeData& genomeData = _genomeData[_currentGeneration][i];
            
            if (!genomeData.isUnderRepair) {
                // Evaluate fitness for non-repair genomes
                FitnessResultType fitness = _fitnessStrategy->evaluate(genome.get_phenotype(), speciationControl);
                
                // Update species assignment
                uint32_t newSpeciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                _genomeData[_currentGeneration][i].speciesId = newSpeciesId;
                
                // Add to fitness results with global index
                _fitnessResults[_currentGeneration].insert({fitness, i});
            } else {
                // For repair genomes, use default fitness and keep existing species
                _fitnessResults[_currentGeneration].insert({FitnessResultType{}, i});
            }
        }

        // Phase 3: Population Analysis - Species grouping and dynamic data update
        auto speciesGrouping = Population::speciesGrouping(_fitnessResults[_currentGeneration], _genomeData[_currentGeneration], _speciesData);

        // Assert that all genomes in species grouping have empty deltas
        for (const auto& [speciesId, indices] : speciesGrouping) {
            for (size_t globalIndex : indices) {
                assert(Operator::emptyDeltas(_population[_currentGeneration][globalIndex]) && 
                       "Species grouping contains genome with uncleared deltas");
            }
        }

        Population::dynamicDataUpdate(_fitnessResults[_currentGeneration], _genomeData[_currentGeneration], _speciesData, _updateParams);
        
        // Phase 4: Elite/Crossover Replacement
        // Plot elites and crossover pairs
        auto eliteIndices = Population::plotElites(speciesGrouping, _eliteParams);
        auto crossoverPairs = Population::plotCrossover(speciesGrouping, _crossoverParams);
        
        // Find eliminated genomes for replacement
        std::vector<size_t> eliminatedIndices;
        for (size_t i = 0; i < _genomeData[_currentGeneration].size(); ++i) {
            if (_genomeData[_currentGeneration][i].isMarkedForElimination) {
                eliminatedIndices.push_back(i);
            }
        }
        
        LOG_DEBUG("Generation {}: Found {} elites, {} crossover pairs, {} eliminated genomes", 
            generation, eliteIndices.size(), crossoverPairs.size(), eliminatedIndices.size());
        
        // Replace eliminated genomes with elites
        size_t replacementIndex = 0;
        for (size_t eliteIndex : eliteIndices) {
            if (replacementIndex < eliminatedIndices.size()) {
                size_t targetIndex = eliminatedIndices[replacementIndex];
                // Copy elite genome
                _population[_currentGeneration][targetIndex] = _population[_currentGeneration][eliteIndex];
                // Reset genome data for elite (not marked for elimination)
                _genomeData[_currentGeneration][targetIndex] = _genomeData[_currentGeneration][eliteIndex];
                _genomeData[_currentGeneration][targetIndex].isMarkedForElimination = false;
                replacementIndex++;
            } else {
                // Add new elite if we have space
                _population[_currentGeneration].push_back(_population[_currentGeneration][eliteIndex]);
                _genomeData[_currentGeneration].push_back(_genomeData[_currentGeneration][eliteIndex]);
                _genomeData[_currentGeneration].back().isMarkedForElimination = false;
            }
        }
        
        // Replace remaining eliminated genomes with crossover offspring
        for (const auto& [parentAIndex, parentBIndex] : crossoverPairs) {
            if (replacementIndex < eliminatedIndices.size()) {
                size_t targetIndex = eliminatedIndices[replacementIndex];
                
                // Get fitness values for crossover
                FitnessResultType fitnessA, fitnessB;
                for (const auto& [fitness, globalIndex] : _fitnessResults[_currentGeneration]) {
                    if (globalIndex == parentAIndex) fitnessA = fitness;
                    if (globalIndex == parentBIndex) fitnessB = fitness;
                }
                
                // Perform crossover
                Genome offspring = Operator::crossover(
                    _population[_currentGeneration][parentAIndex], fitnessA,
                    _population[_currentGeneration][parentBIndex], fitnessB,
                    _crossoverOperatorParams
                );
                
                // Check for cycles
                bool hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                
                if (!hasCycles) {
                    Operator::phenotypeConstruct(offspring);
                }
                
                // Create genome data for crossover offspring
                Population::DynamicGenomeData crossoverData;
                crossoverData.speciesId = _genomeData[_currentGeneration][parentAIndex].speciesId; // Inherit from first parent
                crossoverData.protectionCounter = 0; // Fresh start for crossover
                crossoverData.isUnderRepair = hasCycles;
                crossoverData.isMarkedForElimination = false;
                
                // Replace eliminated genome
                _population[_currentGeneration][targetIndex] = std::move(offspring);
                _genomeData[_currentGeneration][targetIndex] = crossoverData;
                replacementIndex++;
            } else {
                // Add new crossover offspring if we have space
                // (Similar logic as above but push_back instead of replacement)
                break; // For now, just stop if no more elimination slots
            }
        }
        
        // end of generation loop
    }
    
    // Build results directly
    if (_population[_currentGeneration].empty()) {
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
    
    // Extract best genome from fitness results (last in multimap - highest fitness)
    auto bestIt = _fitnessResults[_currentGeneration].rbegin();
    size_t bestIndex = bestIt->second;
    Genome bestGenome = _population[_currentGeneration][bestIndex];
    FitnessResultType bestFitness = bestIt->first;
    
    // Extract final population and fitness values
    std::vector<Genome> finalPopulation = _population[_currentGeneration];
    std::vector<FitnessResultType> finalFitnessValues;
    
    // Build fitness values in the same order as the population
    finalFitnessValues.reserve(_population[_currentGeneration].size());
    for (size_t i = 0; i < _population[_currentGeneration].size(); ++i) {
        // Find fitness for genome at index i
        FitnessResultType fitness{};
        for (const auto& [fitnessResult, globalIndex] : _fitnessResults[_currentGeneration]) {
            if (globalIndex == i) {
                fitness = fitnessResult;
                break;
            }
        }
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
