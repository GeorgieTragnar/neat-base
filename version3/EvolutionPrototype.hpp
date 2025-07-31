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
#include "version3/analysis/SpeciesGrouping.hpp"
#include "version3/population/SpeciesRanking.hpp"
#include "version3/population/SpeciesEquilibriumElimination.hpp"
#include "version3/population/GenomeEquilibriumElimination.hpp"
#include "version3/population/FilterEliminatedIndices.hpp"
#include "version3/population/PlotElites.hpp"
#include "version3/population/PlotCrossover.hpp"
#include "version3/population/Bootstrap.hpp"
#include "version3/population/ResultExtraction.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/population/GenerationTransition.hpp"
#include "version3/phenotype/PhenotypeConstruct.hpp"
#include "version3/phenotype/PhenotypeUpdateWeight.hpp"
#include "version3/phenotype/PhenotypeUpdateNode.hpp"
#include "version3/phenotype/PhenotypeUpdateConnection.hpp"

#include "logger/Logger.hpp"

namespace Evolution {

auto logger = LOGGER("evolution.EvolutionPrototype");

enum class EvolutionOperation {
    ELITE_COPY,
    REPAIR,
    WEIGHT_MUTATION,
    NODE_MUTATION,
    CONNECTION_MUTATION,
    CONNECTION_REACTIVATION
};

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
        uint32_t generationsCompleted
    ) : _bestGenome(std::move(bestGenome)),
        _bestFitness(bestFitness),
        _generationsCompleted(generationsCompleted) {}

    const Genome& getBestGenome() const { return _bestGenome; }
    FitnessResultType getBestFitness() const { return _bestFitness; }
    const std::vector<Genome>& getFinalPopulation() const { return _finalPopulation; }
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
        const Operator::SpeciesEquilibriumParams& speciesEquilibriumParams,
        const Operator::GenomeEquilibriumParams& genomeEquilibriumParams,
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

    EvolutionOperation chooseEvolutionOperation(const Genome& parentGenome, const DynamicGenomeData& parentData);
    void performEvolutionOperation(size_t i);


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
    Operator::SpeciesEquilibriumParams _speciesEquilibriumParams;
    Operator::GenomeEquilibriumParams _genomeEquilibriumParams;
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

    // Crossover pairs storage for delayed execution
    std::vector<std::pair<uint32_t, uint32_t>> _pendingCrossoverPairs;
};

// Template implementation
template<typename FitnessResultType>
EvolutionPrototype<FitnessResultType>::EvolutionPrototype(
    std::shared_ptr<Analysis::FitnessStrategy<FitnessResultType>> fitnessStrategy,
    uint32_t targetPopulationSize,
    const Operator::InitParams& initParams,
    const Operator::SpeciesEquilibriumParams& speciesEquilibriumParams,
    const Operator::GenomeEquilibriumParams& genomeEquilibriumParams,
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
    _speciesEquilibriumParams(speciesEquilibriumParams),
    _genomeEquilibriumParams(genomeEquilibriumParams),
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
    
    // Bootstrap population using the Bootstrap operator
    Operator::bootstrap<FitnessResultType>(
        _populationContainer,
        _globalIndexRegistry,
        _generation,
        _targetPopulationSize,
        [this]() -> Genome {
            return Operator::init(_historyTracker, _initParams);
        },
        [this](const Genome& genome, uint32_t globalIndex) -> DynamicGenomeData {
            // Calculate species assignment and cycle detection
            uint32_t speciesId = Operator::compatibilityDistance(
                genome, 
                _historyTracker, 
                _compatibilityParams
            );
            bool isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
            
            // Create and return metadata
            DynamicGenomeData metadata;
            metadata.speciesId = speciesId;
            metadata.pendingEliminationCounter = 0;
            metadata.isUnderRepair = isUnderRepair;
            metadata.isMarkedForElimination = false;
            metadata.isElite = false;
            metadata.repairAttempts = 0;
            metadata.genomeIndex = globalIndex;
            metadata.parentAIndex = UINT32_MAX;  // Bootstrap genomes have no parents
            metadata.parentBIndex = UINT32_MAX;
            
            return metadata;
        }
    );
}

template<typename FitnessResultType>
EvolutionResults<FitnessResultType> EvolutionPrototype<FitnessResultType>::run(uint32_t maxGenerations) {

    for (; _generation < maxGenerations; ++_generation) {
        
        // TODO: these two should be atomically a single operation
        _populationContainer.clearGenerationFitnessResults(_generation);
        _populationContainer.announceNewGeneration(_generation);

        const size_t populationSize = _populationContainer.getPopulationSize(_generation - 1);
        
        LOG_DEBUG("Generation {}: Starting with {} genomes", _generation, populationSize);

        // Phase 1: 1:1 Evolution Loop (mutations, elites, repairs)
        for (size_t i = 0; i < populationSize; ++i) {
            switch (_globalIndexRegistry.getState(i)) {
                case GenomeState::HotElimination:
                case GenomeState::ColdElimination:
                    Operator::generationTransition(_globalIndexRegistry, i, Operator::GenerationTransitionParams{});
                case GenomeState::ReadyForReplacement:
                    continue;
                case GenomeState::Active:
                    performEvolutionOperation(i);
                    break;
                default:
                    assert(false && "unhandled Genome State");
                    break;
            }
        }
        
        // Phase 2: Execute Crossover (from previous generation's planned pairs)
        if (!_pendingCrossoverPairs.empty()) {
            LOG_DEBUG("Generation {}: Executing {} crossover pairs from previous generation", 
                     _generation, _pendingCrossoverPairs.size());
            
            for (const auto& [parentAIndex, parentBIndex] : _pendingCrossoverPairs) {
                const Genome& parentA = _populationContainer.getGenome(_generation - 1, parentAIndex);
                const Genome& parentB = _populationContainer.getGenome(_generation - 1, parentBIndex);
                const DynamicGenomeData& parentAData = _populationContainer.getGenomeData(_generation - 1, parentAIndex);
                const DynamicGenomeData& parentBData = _populationContainer.getGenomeData(_generation - 1, parentBIndex);
                
                auto [fitnessA, fitnessB] = Operator::fitnessExtraction(
                    _populationContainer, parentAIndex, parentBIndex, _generation - 1);
                
                Operator::genomePlacement(
                    _populationContainer,
                    _globalIndexRegistry,
                    Operator::crossover(parentA, fitnessA, parentB, fitnessB, _crossoverOperatorParams),
                    [&](const Genome& offspring, uint32_t placementIndex) {
                        uint32_t speciesId = Operator::compatibilityDistance(
                            offspring, _historyTracker, _compatibilityParams);
                        bool hasCycles = Operator::hasCycles(offspring, _cycleDetectionParams);
                        DynamicGenomeData metadata = Operator::createCrossoverDynamicData(
                            parentAData, parentBData, parentAIndex, parentBIndex, 
                            speciesId, hasCycles);
                        
                        if (!hasCycles) {
                            Operator::phenotypeConstruct(_populationContainer, placementIndex, _generation);
                        }
                        
                        return metadata;
                    }, _generation, _fitnessStrategy, speciationControl);
            }
            _pendingCrossoverPairs.clear();
        }
        
        // Phase 3: Population Analysis
        auto populationData = Operator::speciesGrouping(_populationContainer, _generation, _speciesData, _globalIndexRegistry);
        
        // Species ranking - calculate performance metrics and assign ordinal rankings
        Operator::speciesRanking(populationData, _populationContainer, _generation, _speciesData);
        
        // Species equilibrium elimination - eliminate worst-performing species
        Operator::speciesEquilibriumElimination(populationData, _speciesData, _populationContainer, _generation, _speciesEquilibriumParams, _globalIndexRegistry);
        
        // Genome equilibrium elimination - eliminate worst-performing genomes within species
        Operator::genomeEquilibriumElimination(populationData, _speciesData, _populationContainer, _generation, _genomeEquilibriumParams, _globalIndexRegistry);
        
        // Filter out eliminated genomes for elite selection and crossover planning
        Operator::filterEliminatedIndices(populationData, _globalIndexRegistry);
        
        // Phase 3.5: Elite Selection - Mark elites for next generation using filtered grouping
        Operator::plotElites(populationData.speciesGrouping, _eliteParams, _populationContainer, _generation, _globalIndexRegistry, _speciesData);
        
        // Phase 4: Crossover Planning - Store pairs for execution in next generation using filtered grouping
        _pendingCrossoverPairs = Operator::plotCrossover(populationData.speciesGrouping, _crossoverParams, _globalIndexRegistry, _populationContainer, _generation, _speciesData);
        
        LOG_DEBUG("Generation {}: Plotted {} crossover pairs for next generation", 
                 _generation, _pendingCrossoverPairs.size());
        
        // end of generation loop
    }
    
    // Use ResultExtraction operator to get best genome index with validation
    uint32_t bestGlobalIndex = Operator::extractBestGenomeIndex(_populationContainer, _generation, _globalIndexRegistry);
    
    if (bestGlobalIndex == UINT32_MAX) {
        // Create a dummy genome for empty population case
        Genome dummyGenome = Operator::init(_historyTracker, _initParams);
        return EvolutionResults<FitnessResultType>(
            std::move(dummyGenome),
            FitnessResultType{},
            maxGenerations
        );
    }
    // Extract best genome and fitness using the validated index
    Genome bestGenome(_populationContainer.getGenome(_generation, bestGlobalIndex));
    
    // Find the fitness for the best genome
    FitnessResultType bestFitness = Operator::fitnessExtraction(_populationContainer, bestGlobalIndex, bestGlobalIndex, _generation + 1).first;
    
    return EvolutionResults<FitnessResultType>(
        std::move(bestGenome),
        bestFitness,
        maxGenerations
    );
}

template<typename FitnessResultType>
void EvolutionPrototype<FitnessResultType>::performEvolutionOperation(size_t i)
{
    const Genome& parentGenome = _populationContainer.getGenome(_generation - 1, i);
    const DynamicGenomeData& parentData = _populationContainer.getGenomeData(_generation - 1, i);

    switch (chooseEvolutionOperation(parentGenome, parentData)) {
        case EvolutionOperation::ELITE_COPY:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Genome(parentGenome), parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.isUnderRepair = false;
                    assert(Operator::genomeEquals(parentGenome, genome) && "Elite genome copy should be identical to parent");
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        case EvolutionOperation::REPAIR:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Genome(parentGenome), parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.isUnderRepair = false;
                    assert(Operator::genomeEquals(parentGenome, genome) && "Elite genome copy should be identical to parent");
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        case EvolutionOperation::WEIGHT_MUTATION:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Operator::weightMutation(parentGenome, _weightMutationParams),
                parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                    data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                    
                    if (!data.isUnderRepair) {
                        Operator::phenotypeUpdateWeight(_populationContainer, i, _generation);
                    }
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        case EvolutionOperation::NODE_MUTATION:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Operator::nodeMutation(parentGenome, _historyTracker, _nodeMutationParams),
                parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                    data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                    
                    if (!data.isUnderRepair) {
                        Operator::phenotypeUpdateNode(_populationContainer, i, _generation);
                    }
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        case EvolutionOperation::CONNECTION_MUTATION:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Operator::connectionMutation(parentGenome, _historyTracker, _connectionMutationParams),
                parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                    data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                    
                    if (!data.isUnderRepair) {
                        Operator::phenotypeUpdateConnection(_populationContainer, i, _generation);
                    }
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        case EvolutionOperation::CONNECTION_REACTIVATION:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Operator::connectionReactivation(parentGenome, _connectionReactivationParams),
                parentData,
                [&](Genome& genome, DynamicGenomeData& data) {
                    data.speciesId = Operator::compatibilityDistance(genome, _historyTracker, _compatibilityParams);
                    data.isUnderRepair = Operator::hasCycles(genome, _cycleDetectionParams);
                    
                    if (!data.isUnderRepair) {
                        Operator::phenotypeUpdateConnection(_populationContainer, i, _generation);
                    }
                }, _fitnessStrategy, speciationControl, _globalIndexRegistry);
            return;
        default: 
            assert(false && "unhandled evolution operation");
            break;
    }
}

template<typename FitnessResultType>
EvolutionOperation EvolutionPrototype<FitnessResultType>::chooseEvolutionOperation(const Genome& parentGenome, const DynamicGenomeData& parentData)
{
    if (parentData.isUnderRepair)
        return EvolutionOperation::REPAIR;
    else if (parentData.isElite)
        return EvolutionOperation::ELITE_COPY;
    
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double random = dist(_rng);
    
    if (random < _mutationParams.weightMutationProbability)
    {
        return EvolutionOperation::WEIGHT_MUTATION;
    } 
    else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability)
    {
        if (Operator::hasActiveConnections(parentGenome))
            return EvolutionOperation::NODE_MUTATION;
    }
    else if (random < _mutationParams.weightMutationProbability + _mutationParams.nodeMutationProbability + _mutationParams.connectionMutationProbability)
    {
        if (Operator::hasPossibleConnections(parentGenome))
            return EvolutionOperation::CONNECTION_MUTATION;
    }
    else
    {
        if (Operator::hasDisabledConnections(parentGenome))
            return EvolutionOperation::CONNECTION_REACTIVATION;
    }
    return EvolutionOperation::WEIGHT_MUTATION;
}

} // namespace Evolution
