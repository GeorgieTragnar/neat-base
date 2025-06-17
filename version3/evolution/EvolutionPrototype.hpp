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
#include "version3/operator/RepairOperator.hpp"
#include "version3/operator/EmptyDeltas.hpp"

// Population management
#include "version3/population/PopulationData.hpp"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/GenerationPlanner.hpp"
#include "version3/population/GenerationPlannerParams.hpp"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/population/ResolveParentIndices.hpp"
#include "version3/population/ReproductiveInstruction.hpp"

#include "../logger/Logger.hpp"

namespace Evolution {

auto logger = LOGGER("evolution.EvolutionPrototype");
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
        const Operator::RepairOperatorParams& repairParams,
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
    Operator::RepairOperatorParams _repairParams;
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
    const Operator::RepairOperatorParams& repairParams,
    uint32_t randomSeed
) : _fitnessStrategy(std::move(fitnessStrategy)),
    _targetPopulationSize(targetPopulationSize),
    _initParams(initParams),
    _plannerParams(plannerParams),
    _updateParams(updateParams),
    _compatibilityParams(compatibilityParams),
    _repairParams(repairParams),
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
    // Find the genome data for the parent
    uint32_t parentPopulationIndex = instruction.globalParentIndices[0];
    const Population::DynamicGenomeData* parentGenomeData = nullptr;
    for (const auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
        if (genomeData.genomeIndex == parentPopulationIndex) {
            parentGenomeData = &genomeData;
            break;
        }
    }
    if (parentGenomeData == nullptr) {
        // static auto logger = LOGGER("evolution.EvolutionPrototype");
        LOG_ERROR("Could not find parent genome data! Looking for parentPopulationIndex={}", parentPopulationIndex);
        LOG_ERROR("Available genomeIndex values in _genomeData[_lastGeneration]:");
        for (const auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
            LOG_ERROR("  genomeIndex: {}", genomeData.genomeIndex);
        }
    }
    assert(parentGenomeData != nullptr && "Could not find parent genome data");
    
    // Assert that no instruction set should have a parent that is under repair
    assert(!parentGenomeData->isUnderRepair && "No instruction set should have a parent under repair - species grouping should exclude such parents");
    
    // Assert that parent genomes for normal operations have empty deltas
    assert(Operator::emptyDeltas(_population[_lastGeneration][parentPopulationIndex]) && "Parent genome 0 has uncleared deltas for normal operation");
    if (instruction.globalParentIndices.size() > 1) {
        assert(Operator::emptyDeltas(_population[_lastGeneration][instruction.globalParentIndices[1]]) && "Parent genome 1 has uncleared deltas for normal operation");
    }
    
    switch (instruction.operationType) {
        case Population::OperationType::PRESERVE: {
            return _population[_lastGeneration][parentPopulationIndex];
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
                    Genome offspring = Operator::weightMutation(_population[_lastGeneration][parentPopulationIndex]
                        , weightParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        // static auto logger = LOGGER("evolution.EvolutionPrototype");
                        LOG_DEBUG("Calling phenotypeUpdateWeight after weight mutation");
                        Operator::phenotypeUpdateWeight(offspring);
                        LOG_DEBUG("phenotypeUpdateWeight completed");
                        // Assert deltas are cleared after phenotype update
                        assert(Operator::emptyDeltas(offspring) && "Weight mutation offspring has uncleared deltas after phenotype update");
                    }
                    return std::move(offspring);
                }
                
                case 1: { // Node mutation
                    NodeGeneAttributes nodeAttribs{ActivationType::SIGMOID}; // TODO: Configure
                    Operator::NodeMutationParams nodeParams(nodeAttribs);
                    Genome offspring = Operator::nodeMutation(_population[_lastGeneration][parentPopulationIndex]
                        , _historyTracker, nodeParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        // static auto logger = LOGGER("evolution.EvolutionPrototype");
                        LOG_DEBUG("Calling phenotypeUpdateNode after node mutation");
                        Operator::phenotypeUpdateNode(offspring);
                        LOG_DEBUG("phenotypeUpdateNode completed");
                        // Assert deltas are cleared after phenotype update
                        assert(Operator::emptyDeltas(offspring) && "Node mutation offspring has uncleared deltas after phenotype update");
                    }
                    return std::move(offspring);
                }
                
                case 2: { // Connection mutation
                    Operator::ConnectionMutationParams connParams(
                        instruction.evolutionParams.mutationRate, 2.0,
                        Operator::ConnectionMutationParams::NetworkTopology::FEED_FORWARD
                    );
                    Genome offspring = Operator::connectionMutation(_population[_lastGeneration][parentPopulationIndex]
                        , _historyTracker, connParams);
                    
                    Operator::CycleDetectionParams cycleParams;
                    hasCycles = Operator::hasCycles(offspring, cycleParams);
                    
                    if (!hasCycles) {
                        // static auto logger = LOGGER("evolution.EvolutionPrototype");
                        LOG_DEBUG("Calling phenotypeUpdateConnection after connection mutation");
                        Operator::phenotypeUpdateConnection(offspring);
                        LOG_DEBUG("phenotypeUpdateConnection completed");
                        // Assert deltas are cleared after phenotype update
                        assert(Operator::emptyDeltas(offspring) && "Connection mutation offspring has uncleared deltas after phenotype update");
                    }
                    return std::move(offspring);
                }
                
                case 3: { // Connection reactivation
                    // TODO: currently disabled
                    // TODO: Check if disabled connections exist first
                    Operator::ConnectionReactivationParams reactParams(
                        Operator::ConnectionReactivationParams::SelectionStrategy::RANDOM
                    );
                    Genome offspring = Operator::connectionReactivation(_population[_lastGeneration][parentPopulationIndex]
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
            // Check if both parents are the same - if so, just copy the first parent
            if (parentPopulationIndex == instruction.globalParentIndices[1]) {
                // static auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Crossover with identical parents ({}), copying first parent instead", parentPopulationIndex);
                return _population[_lastGeneration][parentPopulationIndex];
            }
            
            // Find fitness values for both parents
            const Population::DynamicGenomeData* parent1GenomeData = nullptr;
            FitnessResultType parent0Fitness, parent1Fitness;
            
            for (const auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
                if (genomeData.genomeIndex == parentPopulationIndex) {
                    parent0Fitness = fitness;
                } else if (genomeData.genomeIndex == instruction.globalParentIndices[1]) {
                    parent1Fitness = fitness;
                    parent1GenomeData = &genomeData;
                }
            }
            assert(parent1GenomeData != nullptr && "Could not find parent 1 genome data");
            
            Operator::CrossoverParams crossoverParams;
            Genome offspring = Operator::crossover(
                _population[_lastGeneration][parentPopulationIndex],
                parent0Fitness,
                _population[_lastGeneration][instruction.globalParentIndices[1]],
                parent1Fitness,
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

        // auto logger = LOGGER("evolution.EvolutionPrototype");
        LOG_INFO("generation: {}", generation);

        // Assert that all genomes in last generation with uncleared deltas are properly marked
        if (generation > 0) { // Skip first generation
            for (const auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
                const Genome& genome = _population[_lastGeneration][genomeData.genomeIndex];
                if (!Operator::emptyDeltas(genome)) {
                    assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                           "Genome in last generation has uncleared deltas but is not marked for repair or elimination");
                }
            }
        }

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
                // auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Species {} has {} instruction sets", speciesId, instructionIt->second.size());
                totalGenomesThisGeneration += static_cast<uint32_t>(instructionIt->second.size());
            } else {
                // Species not in instruction sets = new/rediscovered, carry forward
                // auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Species {} has no instructions, carrying forward {} genomes", speciesId, speciesData.currentPopulationSize);
                totalGenomesThisGeneration += speciesData.currentPopulationSize;
            }
        }
        
        // auto logger = LOGGER("evolution.EvolutionPrototype");
        LOG_DEBUG("Generation {}: Calculated total genomes = {}, Current population size = {}", 
            generation, totalGenomesThisGeneration, _population[_currentGeneration].size());
        
        _population[_currentGeneration].clear();
        _population[_currentGeneration].reserve(totalGenomesThisGeneration);
        _genomeData[_currentGeneration].clear();
        
        LOG_DEBUG("Generation {}: Reserved {} slots, current capacity = {}", 
            generation, totalGenomesThisGeneration, _population[_currentGeneration].capacity());

        // Build index-to-iterator mapping once per generation
        std::vector<decltype(_genomeData[_lastGeneration].begin())> indexToIterator;
        indexToIterator.reserve(_genomeData[_lastGeneration].size());
        for (auto it = _genomeData[_lastGeneration].begin(); it != _genomeData[_lastGeneration].end(); ++it) {
            indexToIterator.push_back(it);
        }
        
        // Execute instruction sets and create new generation
        for (const auto& [speciesId, instructions] : _instructionSets[_lastGeneration]) {
            for (const auto& instruction : instructions) {
                // Debug: Log instruction details before createOffspring
                // static auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Executing instruction for species {}: operationType={}, globalParentIndices size={}", 
                    speciesId, static_cast<int>(instruction.operationType), instruction.globalParentIndices.size());
                if (!instruction.globalParentIndices.empty()) {
                    LOG_DEBUG("  Parent 0 index: {}", instruction.globalParentIndices[0]);
                    if (instruction.globalParentIndices.size() > 1) {
                        LOG_DEBUG("  Parent 1 index: {}", instruction.globalParentIndices[1]);
                    }
                }
                
                // Execute the instruction to create offspring
                bool hasCycles = false;
                Genome offspring = createOffspring(instruction, indexToIterator, hasCycles);
                
                // Assert createOffspring contract: either empty deltas or will be marked appropriately
                if (!Operator::emptyDeltas(offspring)) {
                    // If offspring has uncleared deltas, it must be due to cycles (hasCycles=true)
                    assert(hasCycles && "createOffspring returned genome with uncleared deltas but hasCycles=false");
                }

                // Create dynamic genome data for offspring (based on parent)
                // Find the parent's genome data
                Population::DynamicGenomeData genomeData;
                bool foundParentData = false;
                for (const auto& [fitness, parentGenomeData] : _genomeData[_lastGeneration]) {
                    if (parentGenomeData.genomeIndex == instruction.globalParentIndices[0]) {
                        genomeData = parentGenomeData;
                        foundParentData = true;
                        break;
                    }
                }
                assert(foundParentData && "Could not find parent genome data for offspring");
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
                LOG_DEBUG("Adding offspring: population size before = {}, capacity = {}", 
                    _population[_currentGeneration].size(), _population[_currentGeneration].capacity());
                _population[_currentGeneration].push_back(std::move(offspring));
                LOG_DEBUG("Adding offspring: population size after = {}", _population[_currentGeneration].size());
                
                // Assert that any genome with uncleared deltas is properly marked
                const Genome& addedGenome = _population[_currentGeneration].back();
                if (!Operator::emptyDeltas(addedGenome)) {
                    assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                           "Added genome has uncleared deltas but is not marked for repair or elimination");
                }
            }
            for (const auto& [fitness, genomeData] : _genomeData[_currentGeneration]) {
                const Genome& genome = _population[_currentGeneration][genomeData.genomeIndex];
                if (!Operator::emptyDeltas(genome)) {
                    assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                        "Genome entering species grouping has uncleared deltas but is not marked for repair or elimination");
                }
            }
            // instruction set end of loop
        }

        for (const auto& [fitness, genomeData] : _genomeData[_currentGeneration]) {
            const Genome& genome = _population[_currentGeneration][genomeData.genomeIndex];
            if (!Operator::emptyDeltas(genome)) {
                assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                       "Genome entering species grouping has uncleared deltas but is not marked for repair or elimination");
            }
        }
        // 5. Repair parents from last generation that are under repair but not marked for elimination
        // This must happen before carry-forward to ensure genomes have empty deltas
        for (auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
            if (genomeData.isUnderRepair && !genomeData.isMarkedForElimination) {
                // static auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Repairing genome at index {} for species {} (repair attempts: {})", 
                    genomeData.genomeIndex, genomeData.speciesId, genomeData.repairAttempts);
                
                Genome& genomeToRepair = _population[_lastGeneration][genomeData.genomeIndex];
                Genome repairedGenome = Operator::repair(genomeToRepair, genomeData, _repairParams);
                
                // Check if repair was successful
                Operator::CycleDetectionParams cycleParams;
                bool stillHasCycles = Operator::hasCycles(repairedGenome, cycleParams);
                
                // Create genome data for current generation regardless of repair success
                Population::DynamicGenomeData newGenomeData = genomeData;
                newGenomeData.genomeIndex = static_cast<uint32_t>(_population[_currentGeneration].size());
                
                if (!stillHasCycles) {
                    // Repair successful - update fitness and mark as no longer under repair
                    const auto& crepairedGenome = repairedGenome;
                    FitnessResultType newFitness = _fitnessStrategy->evaluate(
                        crepairedGenome.get_phenotype(), speciationControl
                    );
                    
                    // Get species assignment (might have changed after repair)
                    uint32_t newSpeciesId = Operator::compatibilityDistance(
                        repairedGenome, _historyTracker, _compatibilityParams
                    );
                    
                    newGenomeData.speciesId = newSpeciesId;
                    newGenomeData.isUnderRepair = false;
                    
                    _genomeData[_currentGeneration].insert({newFitness, newGenomeData});
                    
                    LOG_DEBUG("Repair successful for genome at index {}, species: {}", 
                        newGenomeData.genomeIndex, newSpeciesId);
                } else {
                    LOG_DEBUG("Repair failed for genome at index {}, will retry next generation", newGenomeData.genomeIndex);
                    
                    // Still under repair, insert with default fitness
                    _genomeData[_currentGeneration].insert({FitnessResultType(), newGenomeData});
                }
                
                // Add the repaired/attempted genome to current generation
                _population[_currentGeneration].push_back(std::move(repairedGenome));
            }
        }
        
        // Handle species that exist in _speciesData but not in instruction sets
        // These are new/rediscovered species that should be copied as-is
        for (const auto& [speciesId, speciesData] : _speciesData) {
            if (_instructionSets[_lastGeneration].find(speciesId) == _instructionSets[_lastGeneration].end()) {
                // auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Processing new/rediscovered species {} with currentPopulationSize={}", 
                    speciesId, speciesData.currentPopulationSize);
                // This species has no instruction sets, copy its genomes from last generation
                for (const auto& [fitness, genomeData] : _genomeData[_lastGeneration]) {
                    if (genomeData.speciesId == speciesId && !genomeData.isUnderRepair) {
                        // Only copy genomes that are not under repair (repair logic handles those above)
                        Genome& sourceGenome = _population[_lastGeneration][genomeData.genomeIndex];
                        
                        // Assert that deltas are empty before copying - genomes should not have uncleared deltas
                        assert(Operator::emptyDeltas(sourceGenome) && "Source genome has uncleared deltas during carry-forward");
                        
                        Genome genomeCopy = sourceGenome;
                        
                        // Create new genome data with updated index
                        Population::DynamicGenomeData newGenomeData = genomeData;
                        newGenomeData.genomeIndex = static_cast<uint32_t>(_population[_currentGeneration].size());
                        
                        // Add to current generation
                        // auto logger = LOGGER("evolution.EvolutionPrototype");
                        LOG_DEBUG("Adding copied genome: population size before = {}, capacity = {}", 
                            _population[_currentGeneration].size(), _population[_currentGeneration].capacity());
                        _population[_currentGeneration].push_back(std::move(genomeCopy));
                        LOG_DEBUG("Adding copied genome: population size after = {}", _population[_currentGeneration].size());
                        _genomeData[_currentGeneration].insert({fitness, newGenomeData});
                    }
                }
            }
        }

        // Population operations
        // 1. Generate instruction sets for next generation based on current species data
        LOG_DEBUG("Generation {}: Before generation planner:", generation);
        for (const auto& [speciesId, speciesData] : _speciesData) {
            LOG_DEBUG("  Species {}: currentPopulationSize={}, instructionSetsSize={}", 
                speciesId, speciesData.currentPopulationSize, speciesData.instructionSetsSize);
        }
        
        for (const auto& [fitness, genomeData] : _genomeData[_currentGeneration]) {
            const Genome& genome = _population[_currentGeneration][genomeData.genomeIndex];
            if (!Operator::emptyDeltas(genome)) {
                assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                       "Genome entering species grouping has uncleared deltas but is not marked for repair or elimination");
            }
        }

        _instructionSets[_currentGeneration] = Population::generationPlanner(_speciesData, _plannerParams);
        
        LOG_DEBUG("Generation {}: After generation planner, instruction sets created:", generation);
        for (const auto& [speciesId, instructions] : _instructionSets[_currentGeneration]) {
            LOG_DEBUG("  Species {}: {} instruction sets", speciesId, instructions.size());
        }

        // Assert that all genomes entering species grouping either have empty deltas or are marked appropriately
        for (const auto& [fitness, genomeData] : _genomeData[_currentGeneration]) {
            const Genome& genome = _population[_currentGeneration][genomeData.genomeIndex];
            if (!Operator::emptyDeltas(genome)) {
                assert((genomeData.isUnderRepair || genomeData.isMarkedForElimination) && 
                       "Genome entering species grouping has uncleared deltas but is not marked for repair or elimination");
            }
        }
        
        // 2. Group current generation by species to get index mappings
        auto speciesGrouping = Population::speciesGrouping(_genomeData[_currentGeneration], _speciesData);
        
        // Debug: Log species grouping results
        LOG_DEBUG("Species grouping results:");
        for (const auto& [speciesId, indices] : speciesGrouping) {
            LOG_DEBUG("  Species {}: {} indices", speciesId, indices.size());
            if (indices.size() <= 10) { // Only log if manageable size
                for (size_t i = 0; i < indices.size(); ++i) {
                    LOG_DEBUG("    index[{}] = {}", i, indices[i]);
                }
            } else if (indices.size() > 10) {
                LOG_DEBUG("    First few indices: {}, {}, {} ... (total {})", 
                    indices[0], indices[1], indices[2], indices.size());
            }
        }
        
        // Assert that all genomes in species grouping have empty deltas
        for (const auto& [speciesId, indices] : speciesGrouping) {
            for (size_t globalIndex : indices) {
                assert(Operator::emptyDeltas(_population[_currentGeneration][globalIndex]) && 
                       "Species grouping contains genome with uncleared deltas");
            }
        }

        // 3. Remove instruction sets for species without valid genomes, then resolve parent indices
        // First, remove instruction sets for species that have no valid genomes in current generation
        for (auto it = _instructionSets[_currentGeneration].begin(); it != _instructionSets[_currentGeneration].end();) {
            auto groupingIt = speciesGrouping.find(it->first);
            if (groupingIt == speciesGrouping.end()) {
                // Species has no valid genomes in current generation - remove entire instruction set
                // static auto logger = LOGGER("evolution.EvolutionPrototype");
                LOG_DEBUG("Removing instruction set for species {} (no valid genomes in current generation)", it->first);
                it = _instructionSets[_currentGeneration].erase(it);
            } else {
                ++it;
            }
        }
        
        // Now resolve parent indices for remaining instruction sets
        for (auto& [speciesId, instructions] : _instructionSets[_currentGeneration]) {
            auto groupingIt = speciesGrouping.find(speciesId);
            assert(groupingIt != speciesGrouping.end() && "Species grouping should exist after cleanup");
            Population::resolveParentIndices(instructions, groupingIt->second);
        }

        // 4. Update dynamic data for current generation (finishing operation of evolution loop)
        LOG_DEBUG("Generation {}: Before dynamic data update:", generation);
        for (const auto& [speciesId, speciesData] : _speciesData) {
            LOG_DEBUG("  Species {}: currentPopulationSize={}, instructionSetsSize={}", 
                speciesId, speciesData.currentPopulationSize, speciesData.instructionSetsSize);
        }
        
        Population::dynamicDataUpdate(_genomeData[_currentGeneration], _speciesData, _instructionSets[_currentGeneration], _updateParams);
        
        LOG_DEBUG("Generation {}: After dynamic data update:", generation);
        for (const auto& [speciesId, speciesData] : _speciesData) {
            LOG_DEBUG("  Species {}: currentPopulationSize={}, instructionSetsSize={}", 
                speciesId, speciesData.currentPopulationSize, speciesData.instructionSetsSize);
        }
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