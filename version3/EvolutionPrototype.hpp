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
#include "version3/population/Init.hpp"
#include "version3/evolution/Crossover.hpp"
#include "version3/population/FitnessExtraction.hpp"
#include "version3/population/CreateCrossoverDynamicData.hpp"
#include "version3/population/GenomePlacement.hpp"
#include "version3/evolution/MutationPlacement.hpp"
#include "version3/evolution/WeightMutation.hpp"
#include "version3/evolution/ConnectionMutation.hpp"
#include "version3/evolution/NodeMutation.hpp"
#include "version3/evolution/ConnectionReactivation.hpp"
#include "version3/validation/CycleDetection.hpp"
#include "version3/analysis/CompatibilityDistance.hpp"
#include "version3/evolution/RepairOperator.hpp"
#include "version3/validation/EmptyDeltas.hpp"
#include "version3/validation/HasDisabledConnections.hpp"
#include "version3/validation/HasActiveConnections.hpp"
#include "version3/validation/HasPossibleConnections.hpp"
#include "version3/validation/GenomeEquals.hpp"

// Population management
#include "version3/data/PopulationContainer.hpp"
#include "version3/data/PopulationData.hpp"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/analysis/SpeciesGrouping.hpp"
#include "version3/population/PlotElites.hpp"
#include "version3/evolution/PlotCrossover.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/population/GenerationTransition.hpp"
#include "version3/phenotype/PhenotypeConstruct.hpp"
#include "version3/phenotype/PhenotypeUpdateWeight.hpp"
#include "version3/phenotype/PhenotypeUpdateNode.hpp"
#include "version3/phenotype/PhenotypeUpdateConnection.hpp"

#include "logger/Logger.hpp"

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
        std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
        uint32_t targetPopulationSize,
        const Operator::InitParams& initParams,
        const Operator::DynamicDataUpdateParams& updateParams,
        const Operator::PlotElitesParams& eliteParams,
        const Operator::PlotCrossoverParams& crossoverParams,
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
    GlobalIndexRegistry _globalIndexRegistry;
    PopulationContainer<FitnessResultType> _populationContainer;
    std::unordered_map<uint32_t, DynamicSpeciesData> _speciesData;

    uint32_t _generation;

    // Strategy and parameters
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> _fitnessStrategy;
    uint32_t _targetPopulationSize;
    Operator::InitParams _initParams;
    Operator::DynamicDataUpdateParams _updateParams;
    Operator::PlotElitesParams _eliteParams;
    Operator::PlotCrossoverParams _crossoverParams;
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
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    uint32_t targetPopulationSize,
    const Operator::InitParams& initParams,
    const Operator::DynamicDataUpdateParams& updateParams,
    const Operator::PlotElitesParams& eliteParams,
    const Operator::PlotCrossoverParams& crossoverParams,
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
        std::vector<std::shared_ptr<const Phenotype>> getChampions() const override {
            return {};
        }
        std::shared_ptr<const Phenotype> getBestChampion() const override {
            return nullptr;
        }
        std::shared_ptr<const Phenotype> getRandomChampion() const override {
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
    
    // Bootstrap population using genomePlacement operator
    for (uint32_t i = 0; i < _targetPopulationSize; ++i) {
        // Create new genome for bootstrap
        Genome newGenome = Operator::init(_historyTracker, _initParams);
        
        // Use genomePlacement with lambda for complete setup
        Operator::genomePlacement<FitnessResultType>(
            _populationContainer,
            _globalIndexRegistry,
            std::move(newGenome),
            [this, &speciationControl](const Genome& genome, size_t placedIndex) -> DynamicGenomeData {
                // Construct phenotype for the placed genome
                Operator::phenotypeConstruct(_populationContainer, placedIndex, 0);
                
                // Calculate species assignment and cycle detection
                uint32_t speciesId = Operator::compatibilityDistance(
                    genome, 
                    _historyTracker, 
                    _compatibilityParams
                );
                bool isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                
                // Create and return complete metadata
                DynamicGenomeData metadata;
                metadata.speciesId = speciesId;
                metadata.pendingEliminationCounter = 0;
                metadata.isUnderRepair = isUnderRepair;
                metadata.isMarkedForElimination = false;
                metadata.parentAIndex = UINT32_MAX;  // Bootstrap genomes have no parents
                metadata.parentBIndex = UINT32_MAX;
                
                return metadata;
            },
            0,  // generation 0
            _fitnessStrategy,
            speciationControl
        );
    }
    
    // Get references for validation and elite selection
    const auto& currentGenomes = _populationContainer.getCurrentGenomes(0);
    const auto& currentGenomeData = _populationContainer.getCurrentGenomeData(0);
    
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
    auto initialSpeciesGrouping = Operator::speciesGrouping(currentFitnessResults, currentGenomeData, _speciesData, _globalIndexRegistry);
    Operator::plotElites(initialSpeciesGrouping, _eliteParams, _globalIndexRegistry, _speciesData);
    
    _generation = 0;
}

template<typename FitnessResultType>
EvolutionResults<FitnessResultType> EvolutionPrototype<FitnessResultType>::run(uint32_t maxGenerations) {
    // Create a simple concrete implementation of SpeciationControlUnit
    class SimpleSpeciationControl : public Analysis::SpeciationControlUnit {
    public:
        std::vector<std::shared_ptr<const Phenotype>> getChampions() const override {
            return {};
        }
        std::shared_ptr<const Phenotype> getBestChampion() const override {
            return nullptr;
        }
        std::shared_ptr<const Phenotype> getRandomChampion() const override {
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
                if (_globalIndexRegistry.getState(i) == GenomeState::Elite) {
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
            const DynamicGenomeData& parentData = lastGenomeData[i];
            
            if (_globalIndexRegistry.getState(i) != GenomeState::Active) {
                // Handle state transitions for eliminated genomes
                auto currentState = _globalIndexRegistry.getState(i);
                
                switch (currentState) {
                    case GenomeState::HotElimination:
                        // One generation in HotElimination -> transition to ColdElimination
                        Operator::generationTransition(_globalIndexRegistry, i, Operator::GenerationTransitionParams{});
                        // LOG_TRACE("Genome {} transitioned: HotElimination -> ColdElimination", i);
                        break;
                        
                    case GenomeState::ColdElimination:
                        // One generation in ColdElimination -> ready for replacement
                        Operator::generationTransition(_globalIndexRegistry, i, Operator::GenerationTransitionParams{});
                        // LOG_TRACE("Genome {} transitioned: ColdElimination -> ReadyForReplacement", i);
                        break;
                        
                    case GenomeState::ReadyForReplacement:
                        // Already ready for replacement, no further transition needed
                        // LOG_TRACE("Genome {} remains ReadyForReplacement", i);
                        break;
                        
                    case GenomeState::Active:
                        // Should not reach here due to outer if condition
                        assert(false && "Logic error: Active genome in non-Active branch");
                        break;
                }
                
                // Copy parent genome/data as placeholder for non-active slots
                currentGenomes[i] = parentGenome;
                currentGenomeData[i] = parentData;
                continue;
            }
            
            // CASE 1: Elite Genomes
            if (_globalIndexRegistry.getState(i) == GenomeState::Elite) {
                Operator::mutationPlacement(_populationContainer, i, _generation,
                    Genome(parentGenome), // Copy without mutation
                    parentData,
                    [&](Genome& genome, DynamicGenomeData& data) {
                        // Elite lambda: minimal updates
                        // Species ID preserved from parent
                        data.isUnderRepair = false;
                        // No phenotype update needed - elite genomes unchanged
                        
                        // Note: Elite fitness handling now managed by mutationPlacement operator
                        
                        assert(Operator::genomeEquals(parentGenome, genome) && "Elite genome copy should be identical to parent");
                    },
                    _fitnessStrategy,
                    speciationControl,
                    _globalIndexRegistry
                );
                continue;
            }
            
            // CASE 2: Repair Genomes
            if (parentData.isUnderRepair) {
                DynamicGenomeData repairData = parentData; // For repair operator
                Genome repaired = Operator::repair(parentGenome, repairData, _repairParams, _globalIndexRegistry, i);
                
                Operator::mutationPlacement(_populationContainer, i, _generation,
                    std::move(repaired),
                    parentData,
                    [&](Genome& genome, DynamicGenomeData& data) {
                        // Repair lambda: check if repair successful
                        bool stillHasCycles = Operator::hasCycles(genome, _cycleDetectionParams);
                        data.isUnderRepair = stillHasCycles;
                        
                        if (!stillHasCycles) {
                            // Repair successful - recalculate species and rebuild phenotype
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            Operator::phenotypeConstruct(_populationContainer, i, _generation);
                            
                            // Note: Fitness evaluation for repaired genomes now handled by mutationPlacement operator
                        }
                        // If still has cycles, keep existing species ID, no phenotype update, no fitness evaluation
                    },
                    _fitnessStrategy,
                    speciationControl,
                    _globalIndexRegistry
                );
                continue;
            }
            
            // CASE 3: Normal Mutations
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            double random = dist(_rng);
            
            if (random < _mutationParams.weightMutationProbability) {
                // Weight Mutation
                Genome mutated = Operator::weightMutation(parentGenome, _weightMutationParams);
                
                Operator::mutationPlacement(_populationContainer, i, _generation,
                    std::move(mutated),
                    parentData,
                    [&](Genome& genome, DynamicGenomeData& data) {
                        // Weight mutation lambda
                        data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                        data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                        
                        if (!data.isUnderRepair) {
                            Operator::phenotypeUpdateWeight(_populationContainer, i, _generation);
                            
                            // Note: Fitness evaluation now handled by mutationPlacement operator
                        }
                    },
                    _fitnessStrategy,
                    speciationControl,
                    _globalIndexRegistry
                );
                
            } else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability) {
                // Node Mutation
                if (Operator::hasActiveConnections(parentGenome)) {
                    Genome mutated = Operator::nodeMutation(parentGenome, _historyTracker, _nodeMutationParams);
                    
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        std::move(mutated),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // Node mutation lambda
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                            
                            if (!data.isUnderRepair) {
                                Operator::phenotypeUpdateNode(_populationContainer, i, _generation);
                                
                                // Note: Fitness evaluation now handled by mutationPlacement operator
                            }
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                } else {
                    // No active connections - copy parent as-is
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        Genome(parentGenome),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // No mutation possible lambda - species and repair status unchanged
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                }
                
            } else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability + _mutationParams.connectionMutationProbability) {
                // Connection Mutation
                if (Operator::hasPossibleConnections(parentGenome)) {
                    Genome mutated = Operator::connectionMutation(parentGenome, _historyTracker, _connectionMutationParams);
                    
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        std::move(mutated),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // Connection mutation lambda
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                            
                            if (!data.isUnderRepair) {
                                Operator::phenotypeUpdateConnection(_populationContainer, i, _generation);
                                
                                // Note: Fitness evaluation now handled by mutationPlacement operator
                            }
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                } else {
                    // Fallback to weight mutation when no connections possible
                    Genome mutated = Operator::weightMutation(parentGenome, _weightMutationParams);
                    
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        std::move(mutated),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // Fallback weight mutation lambda
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                            
                            if (!data.isUnderRepair) {
                                Operator::phenotypeUpdateWeight(_populationContainer, i, _generation);
                                
                                // Note: Fitness evaluation now handled by mutationPlacement operator
                            }
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                }
                
            } else {
                // Connection Reactivation
                if (Operator::hasDisabledConnections(parentGenome)) {
                    Genome mutated = Operator::connectionReactivation(parentGenome, _connectionReactivationParams);
                    
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        std::move(mutated),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // Reactivation lambda
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                            
                            if (!data.isUnderRepair) {
                                Operator::phenotypeUpdateConnection(_populationContainer, i, _generation);
                                
                                // Note: Fitness evaluation now handled by mutationPlacement operator
                            }
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                } else {
                    // Fallback to weight mutation
                    Genome mutated = Operator::weightMutation(parentGenome, _weightMutationParams);
                    
                    Operator::mutationPlacement(_populationContainer, i, _generation,
                        std::move(mutated),
                        parentData,
                        [&](Genome& genome, DynamicGenomeData& data) {
                            // Final fallback weight mutation lambda
                            data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                            data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                            
                            if (!data.isUnderRepair) {
                                Operator::phenotypeUpdateWeight(_populationContainer, i, _generation);
                                
                                // Note: Fitness evaluation now handled by mutationPlacement operator
                            }
                        },
                        _fitnessStrategy,
                        speciationControl,
                        _globalIndexRegistry
                    );
                }
            }
        }

        // Phase 2: Population Analysis - Species grouping and dynamic data update
        // Note: Fitness evaluation now happens within each mutationPlacement lambda
        const auto& lastFitnessResults = _populationContainer.getLastFitnessResults(_generation);
        auto& lastGenomeDataForUpdate = _populationContainer.getLastGenomeData(_generation);
        
        auto speciesGrouping = Operator::speciesGrouping(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, _globalIndexRegistry);

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
            if (_globalIndexRegistry.getState(i) == GenomeState::Elite) {
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
        
        Operator::dynamicDataUpdate(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, speciesGrouping, _updateParams, _globalIndexRegistry);
        
        speciesGrouping = Operator::speciesGrouping(lastFitnessResults, lastGenomeDataForUpdate, _speciesData, _globalIndexRegistry);

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
        Operator::plotElites(speciesGrouping, _eliteParams, _globalIndexRegistry, _speciesData);
        
        // // ELITE_TRACK: Log all genomes marked as Elite after elite selection
        // std::vector<uint32_t> currentElites;
        // std::unordered_map<uint32_t, std::vector<uint32_t>> elitesBySpecies;
        // for (size_t i = 0; i < currentGenomeData.size(); ++i) {
        //     if (_globalIndexRegistry.getState(i) == GenomeState::Elite) {
        //         currentElites.push_back(i);
        //         uint32_t speciesId = currentGenomeData[i].speciesId;
        //         elitesBySpecies[speciesId].push_back(i);
                
        //         // Find fitness for this elite
        //         FitnessResultType eliteFitness{};
        //         for (const auto& [fitness, globalIndex] : currentFitnessResults) {
        //             if (globalIndex == i) {
        //                 eliteFitness = fitness;
        //                 break;
        //             }
        //         }
        //         LOG_DEBUG("ELITE_TRACK: Generation {} - Elite genome {} in species {} with fitness {:.3f}", 
        //                  generation, i, speciesId, eliteFitness.getValue());
        //     }
        // }
        // LOG_DEBUG("ELITE_TRACK: Generation {} - Total {} elites across {} species", 
        //          generation, currentElites.size(), elitesBySpecies.size());
        
        // Phase 4: Crossover Planning
        // Plot crossover pairs for later use
        auto crossoverPairs = Operator::plotCrossover(speciesGrouping, _crossoverParams, _globalIndexRegistry, _speciesData);
        
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
                        if ((state == GenomeState::Active || state == GenomeState::Elite) && 
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
                        if ((state == GenomeState::Active || state == GenomeState::Elite) && 
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
                        if (_globalIndexRegistry.getState(globalIndex) == GenomeState::Elite) {
                            hasElite = true;
                            break;
                        }
                    }
                    assert(hasElite && "Every active species must have at least one elite genome");
                }
            }
        }

        // // Count eliminated genomes for statistics (no longer collecting indices)
        // std::unordered_map<uint32_t, size_t> eliminatedBySpecies;
        // size_t eliminatedUnderRepair = 0;
        // size_t totalEliminated = 0;
        
        // for (size_t i = 0; i < currentGenomeData.size(); ++i) {
        //     const auto& genomeData = currentGenomeData[i];
        //     if (_globalIndexRegistry.getState(i) != GenomeState::Active) {
        //         totalEliminated++;
        //         eliminatedBySpecies[genomeData.speciesId]++;
        //         if (genomeData.isUnderRepair) eliminatedUnderRepair++;
                
        //         // Get fitness for detailed logging
        //         FitnessResultType eliminatedFitness{};
        //         for (const auto& [fitness, globalIndex] : currentFitnessResults) {
        //             if (globalIndex == i) {
        //                 eliminatedFitness = fitness;
        //                 break;
        //             }
        //         }
                
        //         // LOG_TRACE("ELIMINATED GENOME {}: species={}, pendingEliminationCounter={}, fitness={:.3f}, underRepair={}", 
        //         //          i, genomeData.speciesId, genomeData.pendingEliminationCounter, eliminatedFitness.getValue(), genomeData.isUnderRepair);
        //     }
        // }
        
        // LOG_DEBUG("Generation {}: Found {} elites, {} crossover pairs, {} eliminated genomes ({}% population)", 
        //     generation, eliteIndices.size(), crossoverPairs.size(), totalEliminated, 
        //     (totalEliminated * 100.0) / currentGenomeData.size());
        
        // Log current generation state after mutation but before elite/crossover
        // logPopulationState("POST-MUTATION - last", generation, _generation - 1);
        // logPopulationState("POST-MUTATION - current", generation, _generation);
        
        // // Log elimination distribution by species
        // if (!eliminatedBySpecies.empty()) {
        //     std::vector<std::pair<uint32_t, size_t>> elimBySpecies(eliminatedBySpecies.begin(), eliminatedBySpecies.end());
        //     std::sort(elimBySpecies.begin(), elimBySpecies.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
            
        //     std::string elimBySpeciesStr = "[";
        //     for (size_t i = 0; i < elimBySpecies.size(); ++i) {
        //         if (i > 0) elimBySpeciesStr += ", ";
        //         elimBySpeciesStr += fmt::format("{}:{}", elimBySpecies[i].first, elimBySpecies[i].second);
        //     }
        //     elimBySpeciesStr += "]";
        //     // LOG_DEBUG("ELIMINATION DISTRIBUTION: {} under repair, by species: {}", 
        //     //           eliminatedUnderRepair, elimBySpeciesStr);
        // }
        
        // Phase 5: Crossover Replacement - Replace eliminated genomes with crossover offspring
        for (const auto& [parentAIndex, parentBIndex] : crossoverPairs) {
            // Get parent genomes and data
            const auto& parentGenomes = _populationContainer.getGenomes(_generation);
            const auto& parentGenomeData = _populationContainer.getGenomeData(_generation);
            const Genome& parentA = parentGenomes[parentAIndex];
            const Genome& parentB = parentGenomes[parentBIndex];
            const DynamicGenomeData& parentAData = parentGenomeData[parentAIndex];
            const DynamicGenomeData& parentBData = parentGenomeData[parentBIndex];
            
            // Extract fitness values
            auto [fitnessA, fitnessB] = Operator::fitnessExtraction(
                _populationContainer, parentAIndex, parentBIndex, _generation);
            
            // Place crossover result with metadata creation and post-placement logic
            Operator::genomePlacement(
                _populationContainer,
                _globalIndexRegistry,
                Operator::crossover(parentA, fitnessA, parentB, fitnessB, _crossoverOperatorParams),
                [&](const Genome& offspring, size_t placementIndex) {
                    // Create metadata
                    uint32_t speciesId = Operator::compatibilityDistance(
                        offspring, _historyTracker, _compatibilityParams);
                    bool hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                    DynamicGenomeData metadata = Operator::createCrossoverDynamicData(
                        parentAData, parentBData, parentAIndex, parentBIndex, 
                        speciesId, hasCycles);
                    
                    // Post-placement logic (will execute after genome is placed)
                    if (!hasCycles) {
                        Operator::phenotypeConstruct(_populationContainer, placementIndex, _generation);
                    }
                    
                    return metadata;
                },
                _generation
            );
        }
        
        // Log final population composition
        const size_t finalPopSize = currentGenomes.size();
        
        // Validate final population state - Active and Elite genomes need to be valid
        size_t activeGenomeCount = 0;
        size_t eliteGenomeCount = 0;
        size_t nonActiveGenomeCount = 0;
        std::unordered_map<uint32_t, size_t> finalSpeciesDistribution;
        
        for (size_t i = 0; i < currentGenomeData.size(); ++i) {
            const auto& genomeData = currentGenomeData[i];
            auto state = _globalIndexRegistry.getState(i);
            
            if (state == GenomeState::Active) {
                // Validate Active genomes have consistent data
                assert(genomeData.speciesId != UINT32_MAX && "Active genome should have valid species ID");
                finalSpeciesDistribution[genomeData.speciesId]++;
                activeGenomeCount++;
            } else if (state == GenomeState::Elite) {
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
