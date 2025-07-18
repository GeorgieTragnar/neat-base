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
#include "version3/operator/HasPossibleConnections.hpp"
#include "version3/operator/GenomeEquals.hpp"

// Population management
#include "version3/population/PopulationContainer.hpp"
#include "version3/population/PopulationData.hpp"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/population/PlotElites.hpp"
#include "version3/population/PlotCrossover.hpp"
#include "version3/population/GlobalIndexRegistry.hpp"

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

    // Core data containers - triple-buffer architecture
    std::shared_ptr<HistoryTracker> _historyTracker;
    Population::GlobalIndexRegistry _globalIndexRegistry;
    Population::PopulationContainer<FitnessResultType> _populationContainer;
    std::unordered_map<uint32_t, Population::DynamicSpeciesData> _speciesData;

    uint32_t _generation;

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

    // Fitness progression tracking
    std::unordered_map<uint32_t, FitnessResultType> _previousBestFitnessBySpecies;

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
    _historyTracker(std::make_shared<HistoryTracker>()),
    _globalIndexRegistry(0),  // Start with size 0, will grow as needed
    _populationContainer(_globalIndexRegistry) {
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
    _populationContainer.reserveCapacity(_targetPopulationSize);
    
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
    genomeData.parentAIndex = UINT32_MAX;  // Initial genome has no parents
    genomeData.parentBIndex = UINT32_MAX;
    
    // Add initial genome to population
    uint32_t initialGenomeIndex = _populationContainer.push_back(0, std::move(initialGenome), std::move(genomeData));

    // Bootstrap population by replicating initial genome to target size
    auto& currentGenomes = _populationContainer.getCurrentGenomes(0);
    auto& currentGenomeData = _populationContainer.getCurrentGenomeData(0);
    auto& currentFitnessResults = _populationContainer.getCurrentFitnessResults(0);
    
    const Genome& initialGenomeRef = currentGenomes[initialGenomeIndex];
    const Population::DynamicGenomeData& initialDataRef = currentGenomeData[initialGenomeIndex];
    
    // Evaluate fitness for initial genome
    FitnessResultType initialFitness = _fitnessStrategy->evaluate(initialGenomeRef.get_phenotype(), speciationControl);
    currentFitnessResults.insert({initialFitness, initialGenomeIndex});
    
    // Create remaining genomes with proper species assignment
    for (uint32_t i = 1; i < _targetPopulationSize; ++i) {
        // Create new genome with randomized weights
        Genome newGenome = Operator::init(_historyTracker, _initParams);
        Operator::phenotypeConstruct(newGenome);
        
        // Determine proper species assignment for this genome
        uint32_t newSpeciesId = Operator::compatibilityDistance(
            newGenome, 
            _historyTracker, 
            _compatibilityParams
        );
        
        // Create genome metadata with correct species assignment
        Population::DynamicGenomeData newGenomeData;
        newGenomeData.speciesId = newSpeciesId;
        newGenomeData.pendingEliminationCounter = 0;
        newGenomeData.isUnderRepair = false;
        newGenomeData.isMarkedForElimination = false;
        newGenomeData.parentAIndex = UINT32_MAX;  // Initial genomes have no parents
        newGenomeData.parentBIndex = UINT32_MAX;
        
        // Add genome to population
        uint32_t newGenomeIndex = _populationContainer.push_back(0, std::move(newGenome), std::move(newGenomeData));
        
        const Genome& newGenomeC = currentGenomes[newGenomeIndex];
        // Evaluate fitness for this genome
        FitnessResultType genomeFitness = _fitnessStrategy->evaluate(newGenomeC.get_phenotype(), speciationControl);
        currentFitnessResults.insert({genomeFitness, newGenomeIndex});
    }
    
    // Bootstrap validation: verify species assignments are correct
    for (size_t i = 0; i < currentGenomes.size(); ++i) {
        uint32_t storedSpeciesId = currentGenomeData[i].speciesId;
        uint32_t calculatedSpeciesId = Operator::compatibilityDistance(
            currentGenomes[i], 
            _historyTracker, 
            _compatibilityParams
        );
        assert(storedSpeciesId == calculatedSpeciesId && 
               "Bootstrap: stored species ID must match calculated species ID");
    }
    
    // Bootstrap elite selection: establish initial elites for each species
    auto initialSpeciesGrouping = Population::speciesGrouping(currentFitnessResults, currentGenomeData, _speciesData, _globalIndexRegistry);
    Population::plotElites(initialSpeciesGrouping, _eliteParams, _globalIndexRegistry, _speciesData);
    
    _generation = 0;
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
        const auto& genomes = _populationContainer.getGenomes(genNumber);
        const auto& genomeData = _populationContainer.getGenomeData(genNumber);
        const auto& fitnessResults = _populationContainer.getFitnessResults(genNumber);
        
        LOG_DEBUG("{} - Generation {} Population State ({} genomes):", phase, generation, genomes.size());
        for (size_t i = 0; i < genomes.size(); ++i) {
            FitnessResultType fitness{};
            // Find fitness for this genome
            for (const auto& [fit, globalIndex] : fitnessResults) {
                if (globalIndex == i) {
                    fitness = fit;
                    break;
                }
            }
            LOG_DEBUG("  Genome[{}]: species={}, fitness={:.3f}, state={}", 
                     i, genomeData[i].speciesId, fitness.getValue(), 
                     static_cast<int>(_globalIndexRegistry.getState(i)));
        }
    };

    for (uint32_t generation = 0; generation < maxGenerations; ++generation) {
        _generation++;

        // Clear current generation fitness results only
        _populationContainer.clearGenerationFitnessResults(_generation);

        const size_t populationSize = _populationContainer.getLastGenomes(_generation).size();
        
        // Reserve capacity for growth
        if (_populationContainer.getGenerationCapacity(_generation) < populationSize * 2) {
            _populationContainer.reserveCapacity(populationSize * 2);
        }

        LOG_DEBUG("Generation {}: Starting with {} genomes", generation, populationSize);
        
        // Phase 1: 1:1 Evolution - Copy and mutate each genome from last generation
        auto& currentGenomes = _populationContainer.getCurrentGenomes(_generation);
        auto& currentGenomeData = _populationContainer.getCurrentGenomeData(_generation);
        const auto& lastGenomes = _populationContainer.getLastGenomes(_generation);
        const auto& lastGenomeData = _populationContainer.getLastGenomeData(_generation);
        
        // ELITE_TRACK: Verify which genomes retained Elite status from previous generation
        if (generation > 0) {
            std::vector<uint32_t> retainedElites;
            for (size_t i = 0; i < populationSize; ++i) {
                if (_globalIndexRegistry.getState(i) == Population::GenomeState::Elite) {
                    retainedElites.push_back(i);
                    LOG_DEBUG("ELITE_TRACK: Generation {} START - Genome {} retained Elite status in species {}", 
                             generation, i, lastGenomeData[i].speciesId);
                }
            }
            LOG_DEBUG("ELITE_TRACK: Generation {} START - {} genomes retained Elite status from previous generation", 
                     generation, retainedElites.size());
        }
        
        for (size_t i = 0; i < populationSize; ++i) {
            const Genome& parentGenome = lastGenomes[i];
            const Population::DynamicGenomeData& parentData = lastGenomeData[i];
            
            if (_globalIndexRegistry.getState(i) != Population::GenomeState::Active) {
                // Handle state transitions for eliminated genomes
                auto currentState = _globalIndexRegistry.getState(i);
                
                switch (currentState) {
                    case Population::GenomeState::HotElimination:
                        // One generation in HotElimination -> transition to ColdElimination
                        _globalIndexRegistry.transitionToCold(i);
                        // LOG_TRACE("Genome {} transitioned: HotElimination -> ColdElimination", i);
                        break;
                        
                    case Population::GenomeState::ColdElimination:
                        // One generation in ColdElimination -> ready for replacement
                        _globalIndexRegistry.markReadyForReplacement(i);
                        // LOG_TRACE("Genome {} transitioned: ColdElimination -> ReadyForReplacement", i);
                        break;
                        
                    case Population::GenomeState::ReadyForReplacement:
                        // Already ready for replacement, no further transition needed
                        // LOG_TRACE("Genome {} remains ReadyForReplacement", i);
                        break;
                        
                    case Population::GenomeState::Active:
                        // Should not reach here due to outer if condition
                        assert(false && "Logic error: Active genome in non-Active branch");
                        break;
                }
                
                // Copy parent genome/data as placeholder for non-active slots
                currentGenomes[i] = parentGenome;
                currentGenomeData[i] = parentData;
                continue;
            }
            
            // Create corresponding genome data
            Population::DynamicGenomeData offspringData = parentData;
            offspringData.parentAIndex = i;  // Parent index from last generation
            offspringData.parentBIndex = UINT32_MAX;  // Single parent (not crossover)
            
            bool hasCycles = false;
            
            // Copy parent genome
            Genome offspring = parentGenome;
            
            // Check if parent is elite - if so, copy as-is without mutation
            if (_globalIndexRegistry.getState(i) == Population::GenomeState::Elite) {
                // Elite protection: copy as-is without any mutation
                offspring = parentGenome;
                hasCycles = false; // Elites are assumed to be cycle-free
                offspringData.isUnderRepair = false;
                
                // Assert that elite copy is identical
                assert(Operator::genomeEquals(parentGenome, offspring) && "Elite genome copy should be identical to parent");
            } else if (parentData.isUnderRepair) {
                offspring = Operator::repair(offspring, offspringData, _repairParams, _globalIndexRegistry, i);
                
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
                    // Connection mutation - check if new connections are possible
                    if (Operator::hasPossibleConnections(offspring)) {
                        offspring = Operator::connectionMutation(offspring, _historyTracker, _connectionMutationParams);
                        hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        if (!hasCycles) {
                            Operator::phenotypeUpdateConnection(offspring);
                        }
                    } else {
                        // Fallback to weight mutation when no connections possible
                        offspring = Operator::weightMutation(offspring, _weightMutationParams);
                        hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        if (!hasCycles) {
                            Operator::phenotypeUpdateWeight(offspring);
                        }
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
            currentGenomes[i] = std::move(offspring);
            currentGenomeData[i] = offspringData;
        }

        // Phase 2: Fitness Evaluation - Build fitness results multimap
        auto& currentFitnessResults = _populationContainer.getCurrentFitnessResults(_generation);
        
        for (size_t i = 0; i < currentGenomes.size(); ++i) {
            const Genome& genome = currentGenomes[i];
            const Population::DynamicGenomeData& genomeData = currentGenomeData[i];
            
            if (_globalIndexRegistry.getState(i) == Population::GenomeState::ReadyForReplacement) {
                // Skip ReadyForReplacement genomes entirely - they are recycled slots
                continue;
            } else if (_globalIndexRegistry.getState(i) == Population::GenomeState::Elite) {
                // Elite genomes: preserve species assignment and copy fitness from previous generation
                // No need to recalculate compatibility or fitness since elite genomes are unchanged
                
                // Find fitness for this elite from previous generation
                const auto& lastFitnessResults = _populationContainer.getLastFitnessResults(_generation);
                FitnessResultType eliteFitness{};
                bool fitnessFound = false;
                for (const auto& [fitness, globalIndex] : lastFitnessResults) {
                    if (globalIndex == i) {
                        eliteFitness = fitness;
                        fitnessFound = true;
                        break;
                    }
                }
                assert(fitnessFound && "Elite genome should have fitness from previous generation");
                
                // Keep existing species assignment (no compatibility distance recalculation)
                // currentGenomeData[i].speciesId remains unchanged
                
                // Add to fitness results with previous generation's fitness
                currentFitnessResults.insert({eliteFitness, i});
            } else if (!genomeData.isUnderRepair) {
                // Evaluate fitness for non-repair, non-elite genomes
                FitnessResultType fitness = _fitnessStrategy->evaluate(genome.get_phenotype(), speciationControl);
                // Update species assignment
                uint32_t newSpeciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                currentGenomeData[i].speciesId = newSpeciesId;
                
                // Add to fitness results with global index
                currentFitnessResults.insert({fitness, i});
            } else {
                // For repair genomes, use default fitness and keep existing species
                currentFitnessResults.insert({FitnessResultType{}, i});
            }
        }

        // Phase 3: Population Analysis - Species grouping and dynamic data update
        const auto& lastFitnessResults = _populationContainer.getLastFitnessResults(_generation);
        auto& lastGenomeDataForUpdate = _populationContainer.getLastGenomeData(_generation);
        
        auto speciesGrouping = Population::speciesGrouping(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, _globalIndexRegistry);

        // Assert that all genomes in species grouping have empty deltas
        for (const auto& [speciesId, indices] : speciesGrouping) {
            for (size_t globalIndex : indices) {
                assert(Operator::emptyDeltas(lastGenomes[globalIndex]) && 
                       "Species grouping contains genome with uncleared deltas");
            }
        }

        // ELITE_TRACK: Check if any species with Elite genomes are about to be eliminated
        std::unordered_set<uint32_t> speciesWithElites;
        for (size_t i = 0; i < lastGenomeDataForUpdate.size(); ++i) {
            if (_globalIndexRegistry.getState(i) == Population::GenomeState::Elite) {
                speciesWithElites.insert(lastGenomeDataForUpdate[i].speciesId);
            }
        }
        
        std::unordered_set<uint32_t> eliminatedSpeciesIds;
        for (const auto& [speciesId, data] : _speciesData) {
            if (data.isMarkedForElimination) {
                eliminatedSpeciesIds.insert(speciesId);
                if (speciesWithElites.find(speciesId) != speciesWithElites.end()) {
                    LOG_DEBUG("ELITE_TRACK: WARNING - Species {} marked for elimination but contains Elite genomes!", speciesId);
                }
            }
        }
        
        Population::dynamicDataUpdate(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, speciesGrouping, _updateParams, _globalIndexRegistry);
        
        speciesGrouping = Population::speciesGrouping(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, _globalIndexRegistry);

        // DEBUG: Verify all active species appear in the grouping
        for (const auto& [speciesId, speciesData] : _speciesData) {
            if (!speciesData.isMarkedForElimination && speciesData.currentPopulationSize > 0) {
                assert(speciesGrouping.find(speciesId) != speciesGrouping.end() && 
                       "Active species must appear in species grouping after creation");
            }
        }


        // ELITE_TRACK: Verify species elimination didn't affect Elite species
        for (const auto& [speciesId, data] : _speciesData) {
            if (data.isMarkedForElimination && eliminatedSpeciesIds.find(speciesId) == eliminatedSpeciesIds.end()) {
                // This species was newly marked for elimination
                if (speciesWithElites.find(speciesId) != speciesWithElites.end()) {
                    LOG_DEBUG("ELITE_TRACK: CRITICAL - Species {} with Elite genomes was eliminated by dynamicDataUpdate!", speciesId);
                }
            }
        }
        
        // Phase 3.5: Elite Selection - Mark elites in global registry for protection during 1:1 evolution
        Population::plotElites(speciesGrouping, _eliteParams, _globalIndexRegistry, _speciesData);
        
        // ELITE_TRACK: Log all genomes marked as Elite after elite selection
        std::vector<uint32_t> currentElites;
        std::unordered_map<uint32_t, std::vector<uint32_t>> elitesBySpecies;
        for (size_t i = 0; i < currentGenomeData.size(); ++i) {
            if (_globalIndexRegistry.getState(i) == Population::GenomeState::Elite) {
                currentElites.push_back(i);
                uint32_t speciesId = currentGenomeData[i].speciesId;
                elitesBySpecies[speciesId].push_back(i);
                
                // Find fitness for this elite
                FitnessResultType eliteFitness{};
                for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                    if (globalIndex == i) {
                        eliteFitness = fitness;
                        break;
                    }
                }
                LOG_DEBUG("ELITE_TRACK: Generation {} - Elite genome {} in species {} with fitness {:.3f}", 
                         generation, i, speciesId, eliteFitness.getValue());
            }
        }
        LOG_DEBUG("ELITE_TRACK: Generation {} - Total {} elites across {} species", 
                 generation, currentElites.size(), elitesBySpecies.size());
        
        // Phase 4: Crossover Planning
        // Plot crossover pairs for later use
        auto crossoverPairs = Population::plotCrossover(speciesGrouping, _crossoverParams, _globalIndexRegistry, _speciesData);
        
        // NEAT Algorithm Validation: Every active species must have at least one elite
        // Count truly active species (exist in species data, not marked for elimination, and have non-zero population)
        size_t activeSpeciesCount = 0;
        for (const auto& [speciesId, speciesData] : _speciesData) {
            if (!speciesData.isMarkedForElimination && speciesData.currentPopulationSize > 0) {
                activeSpeciesCount++;
            }
        }
        
        // DEBUG: Verify each active species has at least one valid genome in fitness results
        for (const auto& [speciesId, speciesData] : _speciesData) {
            if (!speciesData.isMarkedForElimination && speciesData.currentPopulationSize > 0) {
                // Check if this species has any valid genomes in fitness results
                bool hasValidGenomes = false;
                for (const auto& [fitness, globalIndex] : lastFitnessResults) {
                    if (lastGenomeDataForUpdate[globalIndex].speciesId == speciesId) {
                        auto state = _globalIndexRegistry.getState(globalIndex);
                        if ((state == Population::GenomeState::Active || state == Population::GenomeState::Elite) && 
                            !lastGenomeDataForUpdate[globalIndex].isUnderRepair) {
                            hasValidGenomes = true;
                            break;
                        }
                    }
                }
                assert(hasValidGenomes && "Active species must have at least one valid genome in fitness results");
            }
        }
        
        // DEBUG: Log active species vs grouping details
        LOG_DEBUG("Active species count: {}", activeSpeciesCount);
        LOG_DEBUG("Species grouping size: {}", speciesGrouping.size());
        for (const auto& [speciesId, speciesData] : _speciesData) {
            if (!speciesData.isMarkedForElimination && speciesData.currentPopulationSize > 0) {
                LOG_DEBUG("Active species {}: in grouping = {}, population = {}", speciesId, 
                         speciesGrouping.find(speciesId) != speciesGrouping.end(),
                         speciesData.currentPopulationSize);
            }
        }
        
        // Assert species fitness progression: non-eliminated species should maintain or improve
        for (const auto& [speciesId, speciesData] : _speciesData) {
            // Skip eliminated species and newly discovered species
            if (!speciesData.isMarkedForElimination && 
                _previousBestFitnessBySpecies.find(speciesId) != _previousBestFitnessBySpecies.end()) {
                
                // Find best fitness for this species in current generation
                FitnessResultType currentBest = std::numeric_limits<FitnessResultType>::lowest();
                bool foundGenome = false;
                
                for (const auto& [fitness, globalIndex] : lastFitnessResults) {
                    if (lastGenomeDataForUpdate[globalIndex].speciesId == speciesId) {
                        auto state = _globalIndexRegistry.getState(globalIndex);
                        if ((state == Population::GenomeState::Active || state == Population::GenomeState::Elite) && 
                            !lastGenomeDataForUpdate[globalIndex].isUnderRepair) {
                            if (!foundGenome || fitness.isBetterThan(currentBest)) {
                                currentBest = fitness;
                                foundGenome = true;
                            }
                        }
                    }
                }
                
                if (foundGenome) {
                    FitnessResultType previousBest = _previousBestFitnessBySpecies[speciesId];
                    assert((currentBest.isBetterThan(previousBest) || currentBest.isEqualTo(previousBest)) && 
                           "Non-eliminated species must maintain or improve best fitness");
                }
            }
        }
        
        // Only validate if we have species data from previous generations
        if (activeSpeciesCount > 0) {
            // Validate species grouping completeness against active species count
            assert(speciesGrouping.size() == activeSpeciesCount && 
                   "Every active species should have at least one valid genome for elite selection");
            
            // Per-species validation: every active species must appear in grouping with valid genomes AND elites
            for (const auto& [speciesId, speciesData] : _speciesData) {
                if (!speciesData.isMarkedForElimination && speciesData.currentPopulationSize > 0) {
                    assert(speciesGrouping.find(speciesId) != speciesGrouping.end() && 
                           "Active species must appear in species grouping");
                    assert(!speciesGrouping.at(speciesId).empty() && 
                           "Active species must have at least one valid genome");
                           
                    // Check that at least one genome in this species is Elite
                    bool hasElite = false;
                    for (size_t globalIndex : speciesGrouping.at(speciesId)) {
                        if (_globalIndexRegistry.getState(globalIndex) == Population::GenomeState::Elite) {
                            hasElite = true;
                            break;
                        }
                    }
                    assert(hasElite && "Every active species must have at least one elite genome");
                }
            }
        }

        // Count eliminated genomes for statistics (no longer collecting indices)
        std::unordered_map<uint32_t, size_t> eliminatedBySpecies;
        size_t eliminatedUnderRepair = 0;
        size_t totalEliminated = 0;
        
        for (size_t i = 0; i < currentGenomeData.size(); ++i) {
            const auto& genomeData = currentGenomeData[i];
            if (_globalIndexRegistry.getState(i) != Population::GenomeState::Active) {
                totalEliminated++;
                eliminatedBySpecies[genomeData.speciesId]++;
                if (genomeData.isUnderRepair) eliminatedUnderRepair++;
                
                // Get fitness for detailed logging
                FitnessResultType eliminatedFitness{};
                for (const auto& [fitness, globalIndex] : currentFitnessResults) {
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
        //     generation, eliteIndices.size(), crossoverPairs.size(), totalEliminated, 
        //     (totalEliminated * 100.0) / currentGenomeData.size());
        
        // Log current generation state after mutation but before elite/crossover
        // logPopulationState("POST-MUTATION - last", generation, _generation - 1);
        // logPopulationState("POST-MUTATION - current", generation, _generation);
        
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
        
        // Phase 5: Crossover Replacement - Replace eliminated genomes with crossover offspring
        size_t crossoverReplacements = 0;
        size_t crossoverAdditions = 0;
        size_t crossoverWithCycles = 0;
        
        // LOG_DEBUG("CROSSOVER REPLACEMENT: Starting crossover replacement with {} crossover pairs", 
        //           crossoverPairs.size());
        
        for (const auto& [parentAIndex, parentBIndex] : crossoverPairs) {
            // Get a free index from registry
            uint32_t targetIndex = _globalIndexRegistry.getFreeIndex();
            
            // Get fitness values for crossover
            FitnessResultType fitnessA{}, fitnessB{};
            for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                if (globalIndex == parentAIndex) fitnessA = fitness;
                if (globalIndex == parentBIndex) fitnessB = fitness;
            }
            
            // LOG_TRACE("CROSSOVER REPLACEMENT: Parents {} (species={}, fitness={:.3f}) x {} (species={}, fitness={:.3f}) -> target {}", 
            //          parentAIndex, currentGenomeData[parentAIndex].speciesId, fitnessA.getValue(),
            //          parentBIndex, currentGenomeData[parentBIndex].speciesId, fitnessB.getValue(),
            //          targetIndex);
            
            // Validate crossover parents are not eliminated
            auto parentAState = _globalIndexRegistry.getState(parentAIndex);
            auto parentBState = _globalIndexRegistry.getState(parentBIndex);
            assert((parentAState == Population::GenomeState::Active || parentAState == Population::GenomeState::Elite) && 
                   "Crossover parent A should be in Active or Elite state");
            assert((parentBState == Population::GenomeState::Active || parentBState == Population::GenomeState::Elite) && 
                   "Crossover parent B should be in Active or Elite state");
            
            // Perform crossover
            Genome offspring = Operator::crossover(
                currentGenomes[parentAIndex], fitnessA,
                currentGenomes[parentBIndex], fitnessB,
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
            crossoverData.speciesId = currentGenomeData[parentAIndex].speciesId; // Inherit from first parent
            crossoverData.pendingEliminationCounter = 0; // Fresh start for crossover
            crossoverData.isUnderRepair = hasCycles;
            crossoverData.parentAIndex = parentAIndex;  // First parent from current generation
            crossoverData.parentBIndex = parentBIndex;  // Second parent from current generation
            
            // LOG_TRACE("CROSSOVER OFFSPRING: Created for target {}, species={}, hasCycles={}", 
            //          targetIndex, crossoverData.speciesId, hasCycles);
            
            if (targetIndex == Population::INVALID_INDEX) {
                // No free indices available - add new genome to population
                targetIndex = _populationContainer.push_back(_generation, std::move(offspring), std::move(crossoverData));
                crossoverAdditions++;
            } else {
                // Using existing eliminated slot (already atomically set to Active by getFreeIndex)
                currentGenomes[targetIndex] = std::move(offspring);
                currentGenomeData[targetIndex] = crossoverData;
                crossoverReplacements++;
            }
        }
        
        // LOG_DEBUG("CROSSOVER REPLACEMENT COMPLETE: {} crossover replacements, {} crossover additions, {} with cycles", 
        //           crossoverReplacements, crossoverAdditions, crossoverWithCycles);
        
        // Log final population composition
        const size_t finalPopSize = currentGenomes.size();
        // LOG_DEBUG("GENERATION {} REPLACEMENT SUMMARY: Final pop size: {}, Eliminated: {}, Crossover replacements: {}, Crossover additions: {}", 
        //           generation, finalPopSize, totalEliminated, crossoverReplacements, crossoverAdditions);
        
        // Validate final population state - Active and Elite genomes need to be valid
        size_t activeGenomeCount = 0;
        size_t eliteGenomeCount = 0;
        size_t nonActiveGenomeCount = 0;
        std::unordered_map<uint32_t, size_t> finalSpeciesDistribution;
        
        for (size_t i = 0; i < currentGenomeData.size(); ++i) {
            const auto& genomeData = currentGenomeData[i];
            auto state = _globalIndexRegistry.getState(i);
            
            if (state == Population::GenomeState::Active) {
                // Validate Active genomes have consistent data
                assert(genomeData.speciesId != UINT32_MAX && "Active genome should have valid species ID");
                finalSpeciesDistribution[genomeData.speciesId]++;
                activeGenomeCount++;
            } else if (state == Population::GenomeState::Elite) {
                // Validate Elite genomes have consistent data
                assert(genomeData.speciesId != UINT32_MAX && "Elite genome should have valid species ID");
                finalSpeciesDistribution[genomeData.speciesId]++;
                eliteGenomeCount++;
            } else {
                // Non-Active genomes (HotElimination, ColdElimination, ReadyForReplacement) are acceptable
                nonActiveGenomeCount++;
            }
        }
        
        LOG_DEBUG("GENERATION {} FINAL STATE: {} active genomes, {} elite genomes, {} non-active genomes, {} species active", 
                  generation, activeGenomeCount, eliteGenomeCount, nonActiveGenomeCount, finalSpeciesDistribution.size());
        
        // Log final generation state after all replacements
        // logPopulationState("FINAL-STATE - last", generation, _generation - 1);
        // logPopulationState("FINAL-STATE - current", generation, _generation);
        
        // Update previous generation's best fitness for next iteration
        _previousBestFitnessBySpecies.clear();
        for (const auto& [speciesId, indices] : speciesGrouping) {
            // Find best fitness for this species in current generation
            FitnessResultType bestFitness = std::numeric_limits<FitnessResultType>::lowest();
            bool foundGenome = false;
            
            for (size_t globalIndex : indices) {
                for (const auto& [fitness, fitnessGlobalIndex] : lastFitnessResults) {
                    if (fitnessGlobalIndex == globalIndex) {
                        if (!foundGenome || fitness.isBetterThan(bestFitness)) {
                            bestFitness = fitness;
                            foundGenome = true;
                        }
                        break;
                    }
                }
            }
            
            if (foundGenome) {
                _previousBestFitnessBySpecies[speciesId] = bestFitness;
            }
        }
        
        // end of generation loop
    }
    
    // Build results directly
    const auto& finalGenomes = _populationContainer.getCurrentGenomes(_generation);
    const auto& finalFitnessResults = _populationContainer.getCurrentFitnessResults(_generation);
    
    if (finalGenomes.empty()) {
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
    auto bestIt = finalFitnessResults.rbegin();
    size_t bestIndex = bestIt->second;
    Genome bestGenome = finalGenomes[bestIndex];
    FitnessResultType bestFitness = bestIt->first;
    
    // Extract final population and fitness values
    std::vector<Genome> finalPopulation = finalGenomes;
    std::vector<FitnessResultType> finalFitnessValues;
    
    // Build fitness values in the same order as the population
    finalFitnessValues.reserve(finalGenomes.size());
    for (size_t i = 0; i < finalGenomes.size(); ++i) {
        // Find fitness for genome at index i
        FitnessResultType fitness{};
        for (const auto& [fitnessResult, globalIndex] : finalFitnessResults) {
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
