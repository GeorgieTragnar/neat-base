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

    // Crossover pairs storage for delayed execution
    std::vector<std::pair<size_t, size_t>> _pendingCrossoverPairs;
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
            [this](const Genome& genome, size_t placedIndex) -> DynamicGenomeData {
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
    auto initialSpeciesGrouping = Operator::speciesGrouping(_populationContainer, _generation, _speciesData, _globalIndexRegistry);
    Operator::plotElites(initialSpeciesGrouping, _eliteParams, _populationContainer, _generation, _globalIndexRegistry, _speciesData);
    
    _generation = 0;
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
                    [&](const Genome& offspring, size_t placementIndex) {
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
        auto speciesGrouping = Operator::speciesGrouping(_populationContainer, _generation, _speciesData, _globalIndexRegistry);
        
        // Get fitness results and genome data for dynamic data update
        const auto& fitnessResults = _populationContainer.getCurrentFitnessResults(_generation);
        auto& genomeData = _populationContainer.getCurrentGenomeData(_generation);
        
        Operator::dynamicDataUpdate(fitnessResults, genomeData, _speciesData, speciesGrouping, _updateParams, _globalIndexRegistry);
        
        // Re-run species grouping after dynamic data update (as elimination may have changed)
        speciesGrouping = Operator::speciesGrouping(_populationContainer, _generation, _speciesData, _globalIndexRegistry);
        
        // Phase 3.5: Elite Selection - Mark elites for next generation
        Operator::plotElites(speciesGrouping, _eliteParams, _populationContainer, _generation, _globalIndexRegistry, _speciesData);
        
        // Phase 4: Crossover Planning - Store pairs for execution in next generation
        _pendingCrossoverPairs = Operator::plotCrossover(speciesGrouping, _crossoverParams, _globalIndexRegistry, _populationContainer.getCurrentGenomeData(_generation), _speciesData);
        
        LOG_DEBUG("Generation {}: Plotted {} crossover pairs for next generation", 
                 _generation, _pendingCrossoverPairs.size());
        
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

template<typename FitnessResultType>
void EvolutionPrototype<FitnessResultType>::performEvolutionOperation(size_t i)
{
    const Genome& parentGenome = _populationContainer.getGenome(_generation - 1, i);
    const DynamicGenomeData& parentData = _populationContainer.getGenomeData(_generation - 1, i);

    switch (chooseEvolutionOperation(parentGenome, parentData)) {
        case EvolutionOperation::ELITE_COPY:
            Operator::mutationPlacement(_populationContainer, i, _generation,
                Genome(parentGenome),
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
