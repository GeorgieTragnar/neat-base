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
#include "version3/operator/HasActiveConnections.hpp"

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
    genomeData.pendingEliminationCounter = 0;
    genomeData.isUnderRepair = false;
    genomeData.isMarkedForElimination = false;
    
    // Add to population and genome data (generation 0)
    _population[0].push_back(std::move(initialGenome));
    _genomeData[0].push_back(genomeData);

    const Genome& genome = _population[0][0];
    
    // Evaluate fitness for initial genome
    FitnessResultType initialFitness = _fitnessStrategy->evaluate(genome.get_phenotype(), speciationControl);
    _fitnessResults[0].insert({initialFitness, 0});
    
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

    // Helper function to log complete population state
    auto logPopulationState = [&](const std::string& phase, uint32_t generation, uint32_t genNumber) {
        LOG_DEBUG("{} - Generation {} Population State ({} genomes):", phase, generation, _population[genNumber].size());
        for (size_t i = 0; i < _population[genNumber].size(); ++i) {
            FitnessResultType fitness{};
            // Find fitness for this genome
            for (const auto& [fit, globalIndex] : _fitnessResults[genNumber]) {
                if (globalIndex == i) {
                    fitness = fit;
                    break;
                }
            }
            LOG_DEBUG("  Genome[{}]: species={}, fitness={:.3f}, eliminated={}", 
                     i, _genomeData[genNumber][i].speciesId, fitness.getValue(), 
                     _genomeData[genNumber][i].isMarkedForElimination);
        }
    };

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
                    // Check if genome has any active connections using analytical operator
                    if (Operator::hasActiveConnections(offspring)) {
                        offspring = Operator::nodeMutation(offspring, _historyTracker, _nodeMutationParams);
                        hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        if (!hasCycles) {
                            Operator::phenotypeUpdateNode(offspring);
                        }
                    }
                    // If no active connections, offspring remains unchanged
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
        auto speciesGrouping = Population::speciesGrouping(_fitnessResults[_lastGeneration], _genomeData[_lastGeneration], _speciesData);

        // Assert that all genomes in species grouping have empty deltas
        for (const auto& [speciesId, indices] : speciesGrouping) {
            for (size_t globalIndex : indices) {
                assert(Operator::emptyDeltas(_population[_lastGeneration][globalIndex]) && 
                       "Species grouping contains genome with uncleared deltas");
            }
        }

        Population::dynamicDataUpdate(_fitnessResults[_lastGeneration], _genomeData[_lastGeneration], _speciesData, speciesGrouping, _updateParams);
        
        // Phase 4: Elite/Crossover Replacement
        // Plot elites and crossover pairs
        auto eliteIndices = Population::plotElites(speciesGrouping, _eliteParams);
        auto crossoverPairs = Population::plotCrossover(speciesGrouping, _crossoverParams);
        
        // Find eliminated genomes for replacement
        std::vector<size_t> eliminatedIndices;
        std::unordered_map<uint32_t, size_t> eliminatedBySpecies;
        size_t eliminatedUnderRepair = 0;
        
        for (size_t i = 0; i < _genomeData[_currentGeneration].size(); ++i) {
            const auto& genomeData = _genomeData[_currentGeneration][i];
            if (genomeData.isMarkedForElimination) {
                eliminatedIndices.push_back(i);
                eliminatedBySpecies[genomeData.speciesId]++;
                if (genomeData.isUnderRepair) eliminatedUnderRepair++;
                
                // Get fitness for detailed logging
                FitnessResultType eliminatedFitness{};
                for (const auto& [fitness, globalIndex] : _fitnessResults[_currentGeneration]) {
                    if (globalIndex == i) {
                        eliminatedFitness = fitness;
                        break;
                    }
                }
                
                // LOG_TRACE("ELIMINATED GENOME {}: species={}, pendingEliminationCounter={}, fitness={:.3f}, underRepair={}", 
                //          i, genomeData.speciesId, genomeData.pendingEliminationCounter, eliminatedFitness.getValue(), genomeData.isUnderRepair);
            }
        }
        
        // LOG_DEBUG("Generation {}: Found {} elites, {} crossover pairs, {} eliminated genomes ({}% population)", 
        //     generation, eliteIndices.size(), crossoverPairs.size(), eliminatedIndices.size(), 
        //     (eliminatedIndices.size() * 100.0) / _genomeData[_currentGeneration].size());
        
        // Log current generation state after mutation but before elite/crossover
        // logPopulationState("POST-MUTATION - last", generation, _lastGeneration);
        // logPopulationState("POST-MUTATION - current", generation, _currentGeneration);
        
        // Log elimination distribution by species
        if (!eliminatedBySpecies.empty()) {
            std::vector<std::pair<uint32_t, size_t>> elimBySpecies(eliminatedBySpecies.begin(), eliminatedBySpecies.end());
            std::sort(elimBySpecies.begin(), elimBySpecies.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
            
            std::string elimBySpeciesStr = "[";
            for (size_t i = 0; i < elimBySpecies.size(); ++i) {
                if (i > 0) elimBySpeciesStr += ", ";
                elimBySpeciesStr += fmt::format("{}:{}", elimBySpecies[i].first, elimBySpecies[i].second);
            }
            elimBySpeciesStr += "]";
            // LOG_DEBUG("ELIMINATION DISTRIBUTION: {} under repair, by species: {}", 
            //           eliminatedUnderRepair, elimBySpeciesStr);
        }
        
        // Replace eliminated genomes with elites
        size_t replacementIndex = 0;
        size_t elitesUsedForReplacement = 0;
        size_t elitesAddedToPopulation = 0;
        
        // LOG_DEBUG("ELITE REPLACEMENT: Starting replacement of {} eliminated genomes with {} elites", 
        //           eliminatedIndices.size(), eliteIndices.size());
        
        // Log population state before elite operations
        // logPopulationState("PRE-ELITE-OPS - last", generation, _lastGeneration);
        // logPopulationState("PRE-ELITE-OPS - current", generation, _currentGeneration);
        
        // Debug log: Elite indices and their fitness from last generation
        // LOG_DEBUG("ELITE SELECTION DEBUG: {} elites selected from last generation", eliteIndices.size());
        for (size_t eliteIndex : eliteIndices) {
            FitnessResultType eliteFitness{};
            bool found = false;
            for (const auto& [fitness, globalIndex] : _fitnessResults[_lastGeneration]) {
                if (globalIndex == eliteIndex) {
                    eliteFitness = fitness;
                    found = true;
                    break;
                }
            }
            // LOG_DEBUG("  Elite[{}]: species={}, fitness={:.3f}, found_in_last_gen={}", 
            //          eliteIndex, _genomeData[_lastGeneration][eliteIndex].speciesId, 
            //          eliteFitness.getValue(), found);
        }
        
        for (size_t eliteIndex : eliteIndices) {
            if (replacementIndex < eliminatedIndices.size() && false) {
                size_t targetIndex = eliminatedIndices[replacementIndex];
                
                // Get fitness values for logging
                FitnessResultType eliteFitness{}, eliminatedFitness{};
                for (const auto& [fitness, globalIndex] : _fitnessResults[_lastGeneration]) {
                    if (globalIndex == eliteIndex) eliteFitness = fitness;
                    if (globalIndex == targetIndex) eliminatedFitness = fitness;
                }
                
                
                // Validate replacement is valid
                assert(!_genomeData[_lastGeneration][eliteIndex].isMarkedForElimination && 
                       "Elite genome should not be marked for elimination");
                assert(_genomeData[_currentGeneration][targetIndex].isMarkedForElimination && 
                       "Target genome should be marked for elimination");
                
                // Copy elite genome from last generation
                _population[_currentGeneration][targetIndex] = _population[_lastGeneration][eliteIndex];
                // Reset genome data for elite (copy from last generation, not marked for elimination)
                _genomeData[_currentGeneration][targetIndex] = _genomeData[_lastGeneration][eliteIndex];
                _genomeData[_currentGeneration][targetIndex].isMarkedForElimination = false;
                _fitnessResults[_currentGeneration].insert({eliteFitness, targetIndex});
                replacementIndex++;
                elitesUsedForReplacement++;
            } else {
                // Add new elite if we have space
                FitnessResultType eliteFitness{};
                bool fitnessFound = false;
                // LOG_DEBUG("ELITE FITNESS LOOKUP: Searching for eliteIndex {} in {} fitness entries", 
                        //  eliteIndex, _fitnessResults[_lastGeneration].size());
                for (const auto& [fitness, globalIndex] : _fitnessResults[_lastGeneration]) {
                    // LOG_TRACE("  Checking fitness entry: globalIndex={}, fitness={:.3f}", 
                    //          globalIndex, fitness.getValue());
                    if (globalIndex == eliteIndex) {
                        eliteFitness = fitness;
                        fitnessFound = true;
                        // LOG_DEBUG("  MATCH FOUND: eliteIndex {} has fitness {:.3f}", 
                        //          eliteIndex, fitness.getValue());
                        break;
                    }
                }
                
                if (!fitnessFound) {
                    LOG_ERROR("CRITICAL ERROR: Could not find fitness for eliteIndex {} in _fitnessResults[_lastGeneration]", eliteIndex);
                }
                
                // LOG_TRACE("ELITE ADDITION: Elite {} (species={}, fitness={:.3f}, fitness_found={}) from last generation added to population", 
                //          eliteIndex, _genomeData[_lastGeneration][eliteIndex].speciesId, 
                //          eliteFitness.getValue(), fitnessFound);
                
                _population[_currentGeneration].push_back(_population[_lastGeneration][eliteIndex]);
                _genomeData[_currentGeneration].push_back(_genomeData[_lastGeneration][eliteIndex]);
                _genomeData[_currentGeneration].back().isMarkedForElimination = false;
                _fitnessResults[_currentGeneration].insert({eliteFitness, _population[_currentGeneration].size() - 1});
                elitesAddedToPopulation++;
            }
        }
        
        // LOG_DEBUG("ELITE REPLACEMENT COMPLETE: {} elites used for replacement, {} elites added to population", 
        //           elitesUsedForReplacement, elitesAddedToPopulation);
        
        // Log population state after elite operations
        // logPopulationState("POST-ELITE-OPS - last", generation, _lastGeneration);
        // logPopulationState("POST-ELITE-OPS - current", generation, _currentGeneration);
        
        // Replace remaining eliminated genomes with crossover offspring
        size_t crossoverReplacements = 0;
        size_t crossoverWithCycles = 0;
        
        // LOG_DEBUG("CROSSOVER REPLACEMENT: Starting replacement of {} remaining eliminated genomes with crossover offspring", 
        //           eliminatedIndices.size() - replacementIndex);
        
        for (const auto& [parentAIndex, parentBIndex] : crossoverPairs) {
            if (replacementIndex < eliminatedIndices.size()) {
                size_t targetIndex = eliminatedIndices[replacementIndex];
                
                // Get fitness values for crossover
                FitnessResultType fitnessA{}, fitnessB{}, eliminatedFitness{};
                for (const auto& [fitness, globalIndex] : _fitnessResults[_currentGeneration]) {
                    if (globalIndex == parentAIndex) fitnessA = fitness;
                    if (globalIndex == parentBIndex) fitnessB = fitness;
                    if (globalIndex == targetIndex) eliminatedFitness = fitness;
                }
                
                // LOG_TRACE("CROSSOVER REPLACEMENT: Parents {} (species={}, fitness={:.3f}) x {} (species={}, fitness={:.3f}) -> target {} (fitness={:.3f})", 
                //          parentAIndex, _genomeData[_currentGeneration][parentAIndex].speciesId, fitnessA.getValue(),
                //          parentBIndex, _genomeData[_currentGeneration][parentBIndex].speciesId, fitnessB.getValue(),
                //          targetIndex, eliminatedFitness.getValue());
                
                // Validate crossover parents are not eliminated
                assert(!_genomeData[_currentGeneration][parentAIndex].isMarkedForElimination && 
                       "Crossover parent A should not be marked for elimination");
                assert(!_genomeData[_currentGeneration][parentBIndex].isMarkedForElimination && 
                       "Crossover parent B should not be marked for elimination");
                assert(_genomeData[_currentGeneration][targetIndex].isMarkedForElimination && 
                       "Target genome should be marked for elimination");
                
                // Perform crossover
                Genome offspring = Operator::crossover(
                    _population[_currentGeneration][parentAIndex], fitnessA,
                    _population[_currentGeneration][parentBIndex], fitnessB,
                    _crossoverOperatorParams
                );
                
                // Check for cycles
                bool hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                if (hasCycles) crossoverWithCycles++;
                
                if (!hasCycles) {
                    Operator::phenotypeConstruct(offspring);
                }
                
                // Create genome data for crossover offspring
                Population::DynamicGenomeData crossoverData;
                crossoverData.speciesId = _genomeData[_currentGeneration][parentAIndex].speciesId; // Inherit from first parent
                crossoverData.pendingEliminationCounter = 0; // Fresh start for crossover
                crossoverData.isUnderRepair = hasCycles;
                crossoverData.isMarkedForElimination = false;
                
                // LOG_TRACE("CROSSOVER OFFSPRING: Created for target {}, species={}, hasCycles={}", 
                //          targetIndex, crossoverData.speciesId, hasCycles);
                
                // Replace eliminated genome
                _population[_currentGeneration][targetIndex] = std::move(offspring);
                _genomeData[_currentGeneration][targetIndex] = crossoverData;
                replacementIndex++;
                crossoverReplacements++;
            } else {
                // Add new crossover offspring if we have space
                // (Similar logic as above but push_back instead of replacement)
                // LOG_DEBUG("CROSSOVER EXPANSION: No more elimination slots, stopping crossover expansion");
                break; // For now, just stop if no more elimination slots
            }
        }
        
        // LOG_DEBUG("CROSSOVER REPLACEMENT COMPLETE: {} crossover replacements made, {} with cycles", 
        //           crossoverReplacements, crossoverWithCycles);
        
        // Final replacement validation
        size_t unreplacedEliminations = eliminatedIndices.size() - replacementIndex;
        if (unreplacedEliminations > 0) {
            LOG_ERROR("REPLACEMENT ERROR: {} eliminated genomes were not replaced! Population size will decrease.", 
                     unreplacedEliminations);
            
            // Log details of unreplaced eliminations
            for (size_t i = replacementIndex; i < eliminatedIndices.size(); ++i) {
                size_t unreplacedIndex = eliminatedIndices[i];
                const auto& unreplacedData = _genomeData[_currentGeneration][unreplacedIndex];
                LOG_ERROR("UNREPLACED ELIMINATION {}: species={}, pendingEliminationCounter={}, underRepair={}", 
                         unreplacedIndex, unreplacedData.speciesId, unreplacedData.pendingEliminationCounter, unreplacedData.isUnderRepair);
            }
        }
        
        // Log final population composition
        const size_t finalPopSize = _population[_currentGeneration].size();
        const size_t initialPopSize = finalPopSize - elitesAddedToPopulation; // Approximate initial size
        // LOG_DEBUG("GENERATION {} REPLACEMENT SUMMARY: Pop size {} -> {}, Eliminated: {}, Elite replacements: {}, Crossover replacements: {}, Elite additions: {}", 
        //           generation, initialPopSize, finalPopSize, eliminatedIndices.size(), 
        //           elitesUsedForReplacement, crossoverReplacements, elitesAddedToPopulation);
        
        // Validate final population state
        size_t finalEliminatedCount = 0;
        std::unordered_map<uint32_t, size_t> finalSpeciesDistribution;
        for (size_t i = 0; i < _genomeData[_currentGeneration].size(); ++i) {
            const auto& genomeData = _genomeData[_currentGeneration][i];
            if (genomeData.isMarkedForElimination) {
                finalEliminatedCount++;
                LOG_ERROR("POST-REPLACEMENT ERROR: Genome {} still marked for elimination after replacement phase", i);
            }
            finalSpeciesDistribution[genomeData.speciesId]++;
        }
        
        LOG_DEBUG("GENERATION {} FINAL STATE: {} genomes, {} still marked for elimination, {} species active", 
                  generation, _genomeData[_currentGeneration].size(), finalEliminatedCount, finalSpeciesDistribution.size());
        
        // Log final generation state after all replacements
        // logPopulationState("FINAL-STATE - last", generation, _lastGeneration);
        // logPopulationState("FINAL-STATE - current", generation, _currentGeneration);
        
        // Assert no genomes should be marked for elimination after replacement
        assert(finalEliminatedCount == 0 && "No genomes should remain marked for elimination after replacement");
        
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
