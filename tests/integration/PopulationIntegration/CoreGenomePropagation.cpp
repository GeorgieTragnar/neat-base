#include <gtest/gtest.h>
#include "tests/test_common.h"
#include "tests/test_utilities.h"

#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"
#include "version3/data/PopulationData.hpp"
#include "version3/data/GlobalIndexRegistry.hpp"
#include "version3/data/PopulationContainer.hpp"
#include "version3/operator/CompatibilityDistance.hpp"
#include "version3/operator/Init.hpp"
#include "version3/operator/Crossover.hpp"
#include "version3/operator/WeightMutation.hpp"
#include "version3/operator/NodeMutation.hpp"
#include "version3/operator/ConnectionMutation.hpp"
#include "version3/operator/ConnectionReactivation.hpp"
#include "version3/operator/HasActiveConnections.hpp"
#include "version3/operator/HasPossibleConnections.hpp"
#include "version3/operator/HasDisabledConnections.hpp"
#include "version3/operator/PhenotypeUpdateWeight.hpp"
#include "version3/operator/PhenotypeUpdateNode.hpp"
#include "version3/operator/PhenotypeUpdateConnection.hpp"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/operator/PlotElites.hpp"

#include <memory>
#include <vector>
#include <cmath>
#include <set>
#include <random>

// Equality operators for attribute comparison
bool operator==(const NodeGeneAttributes& a, const NodeGeneAttributes& b) {
    return a.activationType == b.activationType;
}

bool operator==(const ConnectionGeneAttributes& a, const ConnectionGeneAttributes& b) {
    return std::abs(a.weight - b.weight) < 1e-6f && a.enabled == b.enabled;
}

// Core Genome Propagation Integration Test
// Tests the coordination between PopulationContainer, GlobalIndexRegistry,
// HistoryTracker, and CompatibilityDistance without fitness tracking complexity
class CoreGenomePropagationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize core components exactly as specified in integration test plan
        _globalIndexRegistry = std::make_unique<GlobalIndexRegistry>(0);
        _populationContainer = std::make_unique<PopulationContainer<TestFitnessResult>>(*_globalIndexRegistry);
        _historyTracker = std::make_shared<HistoryTracker>();
        _compatibilityParams = std::make_unique<Operator::CompatibilityDistanceParams>(3.0f, 1.0f, 0.4f, 2.0f);
    }
    
    void TearDown() override {
        // Clean up in reverse order
        _compatibilityParams.reset();
        _historyTracker.reset();
        _populationContainer.reset();
        _globalIndexRegistry.reset();
    }

    // Comprehensive genome properties structure for identity tracking
    struct GenomeProperties {
        std::vector<uint32_t> nodeHistoryIDs;
        std::vector<NodeType> nodeTypes;
        std::vector<NodeGeneAttributes> nodeAttributes;
        std::vector<uint32_t> connectionHistoryIDs;
        std::vector<uint32_t> sourceNodeHistoryIDs;
        std::vector<uint32_t> targetNodeHistoryIDs;
        std::vector<ConnectionGeneAttributes> connectionAttributes;
        
        // Byte-for-byte identity comparison
        bool isIdenticalTo(const GenomeProperties& other) const {
            return nodeHistoryIDs == other.nodeHistoryIDs &&
                   nodeTypes == other.nodeTypes &&
                   nodeAttributes == other.nodeAttributes &&
                   connectionHistoryIDs == other.connectionHistoryIDs &&
                   sourceNodeHistoryIDs == other.sourceNodeHistoryIDs &&
                   targetNodeHistoryIDs == other.targetNodeHistoryIDs &&
                   connectionAttributes == other.connectionAttributes;
        }
    };
    
    // Extract all genome properties for baseline comparison
    GenomeProperties extractGenomeProperties(const Genome& genome) {
        GenomeProperties props;
        
        const auto& nodeGenes = genome.get_nodeGenes();
        const auto& connectionGenes = genome.get_connectionGenes();
        
        // Extract node properties
        for (const auto& node : nodeGenes) {
            props.nodeHistoryIDs.push_back(node.get_historyID());
            props.nodeTypes.push_back(node.get_type());
            props.nodeAttributes.push_back(node.get_attributes());
        }
        
        // Extract connection properties
        for (const auto& conn : connectionGenes) {
            props.connectionHistoryIDs.push_back(conn.get_historyID());
            props.sourceNodeHistoryIDs.push_back(nodeGenes[conn.get_sourceNodeIndex()].get_historyID());
            props.targetNodeHistoryIDs.push_back(nodeGenes[conn.get_targetNodeIndex()].get_historyID());
            props.connectionAttributes.push_back(conn.get_attributes());
        }
        
        return props;
    }

    
    // ========================================================================
    // CORE VALIDATION FUNCTIONS
    // ========================================================================
    
    // 1. Global Index Validity
    void validateGlobalIndexConsistency(
        PopulationContainer<TestFitnessResult>& container,
        GlobalIndexRegistry& registry,
        uint32_t gen) {
            
        size_t expectedGenomes = container.getGenerationSize(gen);
        size_t registrySize = registry.getMaxIndex();
        
        // Registry must accommodate all genomes
        EXPECT_GE(registrySize, expectedGenomes) << "Registry too small for generation " << gen;
        
        // All indices should be valid and accessible
        for (size_t i = 0; i < expectedGenomes; ++i) {
            EXPECT_NO_THROW(registry.getState(i)) << "Invalid global index " << i << " in generation " << gen;
        }
        
        // Test buffer rotation and index mapping
        // Verify (generation % 3) mapping works correctly
        uint32_t bufferIndex = gen % 3;
        auto& currentGenomes = container.getGenomes(bufferIndex);
        auto& currentGeneration = container.getCurrentGenomes(gen);
        EXPECT_EQ(&currentGenomes, &currentGeneration) << "Buffer mapping failed for generation " << gen;
    }
    
    // 2. Innovation Number Consistency  
    void validateHistoryTrackerConsistency(
        std::shared_ptr<HistoryTracker> historyTracker,
        uint32_t gen) {
            
        // Test: Same structures get same innovation numbers across generations
        // Verify innovation uniqueness and consistency - same calls should return same IDs
        uint32_t inputNodeId1 = historyTracker->get_input(0);
        uint32_t inputNodeId2 = historyTracker->get_input(0);
        uint32_t outputNodeId1 = historyTracker->get_output(0);
        uint32_t outputNodeId2 = historyTracker->get_output(0);
        uint32_t connectionId1 = historyTracker->get_connection(inputNodeId1, outputNodeId1);
        uint32_t connectionId2 = historyTracker->get_connection(inputNodeId1, outputNodeId1);
        
        // Verify consistency within the same generation
        EXPECT_EQ(inputNodeId1, inputNodeId2) << "Input innovation not consistent in generation " << gen;
        EXPECT_EQ(outputNodeId1, outputNodeId2) << "Output innovation not consistent in generation " << gen;
        EXPECT_EQ(connectionId1, connectionId2) << "Connection innovation not consistent in generation " << gen;
        
        // Verify innovation numbers are reasonable (> 0)
        EXPECT_GT(inputNodeId1, 0u) << "Input node innovation should be positive in generation " << gen;
        EXPECT_GT(outputNodeId1, 0u) << "Output node innovation should be positive in generation " << gen;
        EXPECT_GT(connectionId1, 0u) << "Connection innovation should be positive in generation " << gen;
    }
    
    // 3. Species Assignment Consistency Validation
    void validateSpeciesConsistency(uint32_t gen, const std::vector<std::vector<uint32_t>>& baselineSpeciesIDs, 
                                   size_t numCopiedGenomes) {
        auto& currentGenomes = _populationContainer->getCurrentGenomes(gen);
        auto& currentGenomeData = _populationContainer->getCurrentGenomeData(gen);
        
        // Validate that species assignments are consistent only for copied genomes
        // Mutated genomes are allowed to have different species assignments
        for (size_t i = 0; i < numCopiedGenomes && i < currentGenomes.size() && i < baselineSpeciesIDs[gen].size(); ++i) {
            const auto& genome = currentGenomes[i];
            uint32_t currentSpeciesId = Operator::compatibilityDistance(genome, _historyTracker, *_compatibilityParams);
            uint32_t storedSpeciesId = currentGenomeData[i].speciesId;
            uint32_t baselineSpeciesId = baselineSpeciesIDs[gen][i];
            
            // Verify consistency between calculated, stored, and baseline species IDs for copied genomes
            EXPECT_EQ(currentSpeciesId, baselineSpeciesId) 
                << "Calculated species ID changed for copied genome " << i << " in generation " << gen;
            EXPECT_EQ(storedSpeciesId, baselineSpeciesId) 
                << "Stored species ID inconsistent for copied genome " << i << " in generation " << gen;
        }
        
        // For mutated genomes, just verify that calculated and stored species IDs match
        for (size_t i = numCopiedGenomes; i < currentGenomes.size(); ++i) {
            const auto& genome = currentGenomes[i];
            uint32_t currentSpeciesId = Operator::compatibilityDistance(genome, _historyTracker, *_compatibilityParams);
            uint32_t storedSpeciesId = currentGenomeData[i].speciesId;
            
            EXPECT_EQ(currentSpeciesId, storedSpeciesId) 
                << "Calculated and stored species IDs don't match for mutated genome " << i << " in generation " << gen;
        }
    }
    
    // 4. Comprehensive Genome Consistency Validation
    void validateAllGenomeConsistency(uint32_t gen, const std::vector<std::vector<GenomeProperties>>& baselineProperties,
                                     size_t numCopiedGenomes) {
        auto& currentGenomes = _populationContainer->getCurrentGenomes(gen);
        
        // Validate that copied genomes match their baselines from the previous generation
        if (gen > 0 && numCopiedGenomes > 0) {
            for (size_t i = 0; i < numCopiedGenomes && i < currentGenomes.size(); ++i) {
                const auto& currentGenome = currentGenomes[i];
                GenomeProperties currentProps = extractGenomeProperties(currentGenome);
                
                // Copied genomes should match their baseline from the previous generation
                if (i < baselineProperties[gen - 1].size()) {
                    const auto& baseline = baselineProperties[gen - 1][i];
                    EXPECT_TRUE(currentProps.isIdenticalTo(baseline)) 
                        << "Copied genome " << i << " changed from its baseline in generation " << gen;
                }
            }
        }
        
        // For mutated genomes, we don't validate against baselines since they're expected to change
        // We just verify they have valid structure (done by other validation functions)
    }

private:
    std::unique_ptr<GlobalIndexRegistry> _globalIndexRegistry;
    std::unique_ptr<PopulationContainer<TestFitnessResult>> _populationContainer;
    std::shared_ptr<HistoryTracker> _historyTracker;
    std::unique_ptr<Operator::CompatibilityDistanceParams> _compatibilityParams;
    
    // Helper function to apply mutations with connection validation (similar to EvolutionPrototype)
    struct MutationProbabilities {
        double weightMutation = 0.6;
        double nodeMutation = 0.15;
        double connectionMutation = 0.2;
        double connectionReactivation = 0.05;
    };
    
    Genome applyMutationWithValidation(
        const Genome& parent,
        std::shared_ptr<HistoryTracker> historyTracker,
        const MutationProbabilities& mutationProbs,
        const Operator::WeightMutationParams& weightMutationParams,
        const Operator::NodeMutationParams& nodeMutationParams,
        const Operator::ConnectionMutationParams& connectionMutationParams,
        const Operator::ConnectionReactivationParams& connectionReactivationParams,
        std::mt19937& rng) {
        
        // Copy parent genome
        Genome offspring = parent;
        
        // Generate random number for mutation selection
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double random = dist(rng);
        
        if (random < mutationProbs.weightMutation) {
            // Weight mutation (always possible)
            offspring = Operator::weightMutation(offspring, weightMutationParams);
            Operator::phenotypeUpdateWeight(offspring);
        } else if (random < mutationProbs.weightMutation + mutationProbs.nodeMutation) {
            // Node mutation - check if genome has any active connections
            if (Operator::hasActiveConnections(offspring)) {
                offspring = Operator::nodeMutation(offspring, historyTracker, nodeMutationParams);
                Operator::phenotypeUpdateNode(offspring);
            } else {
                // Fallback to weight mutation if no active connections
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        } else if (random < mutationProbs.weightMutation + mutationProbs.nodeMutation + mutationProbs.connectionMutation) {
            // Connection mutation - check if new connections are possible
            if (Operator::hasPossibleConnections(offspring)) {
                offspring = Operator::connectionMutation(offspring, historyTracker, connectionMutationParams);
                Operator::phenotypeUpdateConnection(offspring);
            } else {
                // Fallback to weight mutation if no connections possible
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        } else {
            // Connection reactivation - check if genome has any disabled connections
            if (Operator::hasDisabledConnections(offspring)) {
                offspring = Operator::connectionReactivation(offspring, connectionReactivationParams);
                Operator::phenotypeUpdateConnection(offspring);
            } else {
                // Fallback to weight mutation if no disabled connections
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        }
        
        return offspring;
    }
};

// ============================================================================
// EVOLUTIONARY FOUNDATION FIXTURE - Provides Evolved State for Higher-Level Tests
// ============================================================================

// Base fixture that runs Level 1 evolution and provides evolved state to derived tests
// This allows Level 2+ tests to start with realistic evolved populations instead of
// duplicating the complex evolutionary setup logic
class EvolutionaryFoundationFixture : public CoreGenomePropagationTest {
protected:
    void SetUpEvolutionaryFoundation() {
        // Run the complete Level 1 evolutionary process to establish realistic foundation
        runCoreGenomePropagationEvolution();
        
        // Capture final evolved state for use by derived test classes
        _foundationGeneration = GENERATIONS - 1;
        _evolvedPopulationSize = _populationContainer->getGenerationSize(_foundationGeneration);
        
        // Analyze the evolved species distribution for Level 2+ testing
        analyzeEvolvedSpeciesDistribution();
        
        std::cout << "\n=== EVOLUTIONARY FOUNDATION ESTABLISHED ===" << std::endl;
        std::cout << "Foundation generation: " << _foundationGeneration << std::endl;
        std::cout << "Population size: " << _evolvedPopulationSize << " genomes" << std::endl;
        std::cout << "Species diversity: " << _evolvedSpeciesMap.size() << " species" << std::endl;
        
        // Log species distribution for diagnostic purposes
        for (const auto& [speciesId, genomeIndices] : _evolvedSpeciesMap) {
            std::cout << "Species " << speciesId << ": " << genomeIndices.size() << " genomes" << std::endl;
        }
        std::cout << "=============================================" << std::endl;
    }

private:
    // Level 1 evolutionary constants (matching the original test)
    static constexpr size_t INITIAL_POPULATION = 7;
    static constexpr size_t GENERATIONS = 8;
    static constexpr size_t NEW_PER_GENERATION = 5;
    
    void runCoreGenomePropagationEvolution() {
        // Exact implementation of Level 1 evolutionary process
        // This is extracted from the original TEST_F to make it reusable
        
        // Create Init operator parameters for consistent genome creation
        std::vector<NodeGeneAttributes> inputAttributes = {{ActivationType::NONE}};
        std::vector<NodeGeneAttributes> outputAttributes = {{ActivationType::SIGMOID}};
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttributes;
        biasAttributes[0] = {1.0f, true};
        
        Operator::InitParams initParams(inputAttributes, outputAttributes, biasAttributes, 
                                      Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
        
        // Store baseline properties for all genomes across generations
        _baselineProperties.resize(GENERATIONS);
        _baselineSpeciesIDs.resize(GENERATIONS);
        
        // Mutation probability parameters (similar to EvolutionPrototype)
        MutationProbabilities mutationProbs;
        mutationProbs.weightMutation = 0.6;
        mutationProbs.nodeMutation = 0.15;
        mutationProbs.connectionMutation = 0.2;
        mutationProbs.connectionReactivation = 0.05;
        
        // Mutation operator parameters
        Operator::WeightMutationParams weightMutationParams(0.8, 0.1, 0.1, 2.0, Operator::WeightMutationParams::MutationType::MIXED);
        Operator::NodeMutationParams nodeMutationParams({ActivationType::TANH});
        Operator::ConnectionMutationParams connectionMutationParams(1.0, 2.0, Operator::ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
        Operator::ConnectionReactivationParams connectionReactivationParams(Operator::ConnectionReactivationParams::SelectionStrategy::RANDOM);
        
        // Random number generator for mutations
        std::random_device rd;
        std::mt19937 gen_rand(rd());
        
        // ========================================================================
        // GENERATION 0 - Create Initial Population Using Init Operator
        // ========================================================================
        
        _baselineProperties[0].reserve(INITIAL_POPULATION);
        _baselineSpeciesIDs[0].reserve(INITIAL_POPULATION);
        
        // Create initial population using Init operator
        for (size_t i = 0; i < INITIAL_POPULATION; ++i) {
            Genome genome = Operator::init(_historyTracker, initParams);
            GenomeProperties props = extractGenomeProperties(genome);
            uint32_t speciesId = Operator::compatibilityDistance(genome, _historyTracker, *_compatibilityParams);
            
            Population::DynamicGenomeData genomeData;
            genomeData.speciesId = speciesId;
            genomeData.pendingEliminationCounter = 0;
            genomeData.repairAttempts = 0;
            genomeData.isUnderRepair = false;
            genomeData.isMarkedForElimination = false;
            genomeData.genomeIndex = UINT32_MAX;
            genomeData.parentAIndex = UINT32_MAX;
            genomeData.parentBIndex = UINT32_MAX;
            
            _baselineProperties[0].push_back(props);
            _baselineSpeciesIDs[0].push_back(speciesId);
            
            uint32_t index = _populationContainer->push_back(0, std::move(genome), genomeData);
            assert(index == i);
        }
        
        // Validate generation 0 setup
        assert(_populationContainer->getGenerationSize(0) == INITIAL_POPULATION);
        validateGlobalIndexConsistency(*_populationContainer, *_globalIndexRegistry, 0);
        validateHistoryTrackerConsistency(_historyTracker, 0);
        
        // ========================================================================
        // EVOLUTIONARY LOOP - Generations 1 through GENERATIONS-1
        // ========================================================================
        
        for (uint32_t gen = 1; gen < GENERATIONS; ++gen) {
            size_t prevGenSize = _populationContainer->getGenerationSize(gen - 1);
            if (prevGenSize == 0) continue;
            
            _baselineProperties[gen].reserve(prevGenSize + NEW_PER_GENERATION);
            _baselineSpeciesIDs[gen].reserve(prevGenSize + NEW_PER_GENERATION);
            
            // Step 1: Copy ALL genomes from last generation to current generation
            auto& lastGenomes = _populationContainer->getLastGenomes(gen);
            auto& lastGenomeData = _populationContainer->getLastGenomeData(gen);
            auto& currentGenomes = _populationContainer->getCurrentGenomes(gen);
            auto& currentGenomeData = _populationContainer->getCurrentGenomeData(gen);
            
            currentGenomes = lastGenomes;
            currentGenomeData = lastGenomeData;
            
            // Store baselines for all copied genomes
            for (size_t i = 0; i < prevGenSize; ++i) {
                if (i < _baselineProperties[gen - 1].size()) {
                    _baselineProperties[gen].push_back(_baselineProperties[gen - 1][i]);
                    _baselineSpeciesIDs[gen].push_back(_baselineSpeciesIDs[gen - 1][i]);
                }
            }
            
            // Step 2: Create new genomes through mutations
            size_t startParentIdx = (prevGenSize >= NEW_PER_GENERATION) ? (prevGenSize - NEW_PER_GENERATION) : 0;
            
            for (size_t newIdx = 0; newIdx < NEW_PER_GENERATION; ++newIdx) {
                size_t parentIdx = startParentIdx + newIdx;
                if (parentIdx >= prevGenSize) {
                    parentIdx = prevGenSize - 1;
                }
                
                const auto& parentGenome = lastGenomes[parentIdx];
                
                // Apply mutations with connection validation
                Genome offspring = applyMutationWithValidation(
                    parentGenome, 
                    _historyTracker, 
                    mutationProbs,
                    weightMutationParams,
                    nodeMutationParams,
                    connectionMutationParams,
                    connectionReactivationParams,
                    gen_rand);
                
                uint32_t newSpeciesId = Operator::compatibilityDistance(offspring, _historyTracker, *_compatibilityParams);
                
                Population::DynamicGenomeData newGenomeData;
                newGenomeData.speciesId = newSpeciesId;
                newGenomeData.pendingEliminationCounter = 0;
                newGenomeData.repairAttempts = 0;
                newGenomeData.isUnderRepair = false;
                newGenomeData.isMarkedForElimination = false;
                newGenomeData.genomeIndex = UINT32_MAX;
                newGenomeData.parentAIndex = parentIdx;
                newGenomeData.parentBIndex = UINT32_MAX;
                
                GenomeProperties newProps = extractGenomeProperties(offspring);
                _baselineProperties[gen].push_back(newProps);
                _baselineSpeciesIDs[gen].push_back(newSpeciesId);
                
                _populationContainer->push_back(gen, std::move(offspring), newGenomeData);
            }
            
            // Validate this generation
            validateGlobalIndexConsistency(*_populationContainer, *_globalIndexRegistry, gen);
            validateHistoryTrackerConsistency(_historyTracker, gen);
            validateSpeciesConsistency(gen, _baselineSpeciesIDs, prevGenSize);
            validateAllGenomeConsistency(gen, _baselineProperties, prevGenSize);
            
            size_t expectedSize = prevGenSize + NEW_PER_GENERATION;
            assert(_populationContainer->getGenerationSize(gen) == expectedSize);
        }
    }
    
    void analyzeEvolvedSpeciesDistribution() {
        auto& finalGenomes = _populationContainer->getCurrentGenomes(_foundationGeneration);
        auto& finalGenomeData = _populationContainer->getCurrentGenomeData(_foundationGeneration);
        
        // Build species distribution map for use by derived test classes
        _evolvedSpeciesMap.clear();
        for (size_t i = 0; i < finalGenomes.size(); ++i) {
            uint32_t speciesId = finalGenomeData[i].speciesId;
            _evolvedSpeciesMap[speciesId].push_back(i);
        }
        
        // Validate that we have species diversity (important for Level 2+ testing)
        assert(!_evolvedSpeciesMap.empty() && "Evolved population should have at least one species");
    }

protected:
    // Evolved state accessible to derived test classes
    uint32_t _foundationGeneration = 0;
    size_t _evolvedPopulationSize = 0;
    std::unordered_map<uint32_t, std::vector<size_t>> _evolvedSpeciesMap;
    
    // Baseline data from evolutionary process (for consistency validation)
    std::vector<std::vector<GenomeProperties>> _baselineProperties;
    std::vector<std::vector<uint32_t>> _baselineSpeciesIDs;
    
    // Helper struct for mutation probabilities (copied from Level 1 test)
    struct MutationProbabilities {
        double weightMutation = 0.6;
        double nodeMutation = 0.15;
        double connectionMutation = 0.2;
        double connectionReactivation = 0.05;
    };
    
    // Helper function for applying mutations (copied from Level 1 test)
    Genome applyMutationWithValidation(
        const Genome& parent,
        std::shared_ptr<HistoryTracker> historyTracker,
        const MutationProbabilities& mutationProbs,
        const Operator::WeightMutationParams& weightMutationParams,
        const Operator::NodeMutationParams& nodeMutationParams,
        const Operator::ConnectionMutationParams& connectionMutationParams,
        const Operator::ConnectionReactivationParams& connectionReactivationParams,
        std::mt19937& rng) {
        
        Genome offspring = parent;
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double random = dist(rng);
        
        if (random < mutationProbs.weightMutation) {
            offspring = Operator::weightMutation(offspring, weightMutationParams);
            Operator::phenotypeUpdateWeight(offspring);
        } else if (random < mutationProbs.weightMutation + mutationProbs.nodeMutation) {
            if (Operator::hasActiveConnections(offspring)) {
                offspring = Operator::nodeMutation(offspring, historyTracker, nodeMutationParams);
                Operator::phenotypeUpdateNode(offspring);
            } else {
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        } else if (random < mutationProbs.weightMutation + mutationProbs.nodeMutation + mutationProbs.connectionMutation) {
            if (Operator::hasPossibleConnections(offspring)) {
                offspring = Operator::connectionMutation(offspring, historyTracker, connectionMutationParams);
                Operator::phenotypeUpdateConnection(offspring);
            } else {
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        } else {
            if (Operator::hasDisabledConnections(offspring)) {
                offspring = Operator::connectionReactivation(offspring, connectionReactivationParams);
                Operator::phenotypeUpdateConnection(offspring);
            } else {
                offspring = Operator::weightMutation(offspring, weightMutationParams);
                Operator::phenotypeUpdateWeight(offspring);
            }
        }
        
        return offspring;
    }
};

// ============================================================================
// FITNESS TRACKING TEST FIXTURE - Level 2 Integration Tests
// ============================================================================

// Level 2 fixture that inherits evolved population from Level 1 and adds fitness tracking
// This tests fitness result storage, elite selection, and elite preservation mechanisms
class FitnessTrackingTest : public EvolutionaryFoundationFixture {
protected:
    void SetUp() override {
        // First establish the base components
        EvolutionaryFoundationFixture::SetUp();
        
        // Then run Level 1 evolution to create realistic foundation
        SetUpEvolutionaryFoundation();
        
        // Initialize fitness tracking components on top of evolved population
        initializeFitnessTrackingComponents();
        
        std::cout << "Level 2 Fitness Tracking Test ready with evolved foundation" << std::endl;
    }

private:
    void initializeFitnessTrackingComponents() {
        // Initialize elite selection parameters
        _eliteParams = std::make_unique<Operator::PlotElitesParams>(0.20, 1, 3); // 20% elites, min 1, max 3
        
        // Initialize species data map (required for elite selection)
        initializeSpeciesDataFromEvolvedPopulation();
        
        std::cout << "Fitness tracking components initialized" << std::endl;
    }
    
    void initializeSpeciesDataFromEvolvedPopulation() {
        // Create species data entries for all evolved species
        for (const auto& [speciesId, genomeIndices] : _evolvedSpeciesMap) {
            DynamicSpeciesData speciesData;
            speciesData.currentPopulationSize = genomeIndices.size();
            speciesData.isMarkedForElimination = false;
            // Initialize other species data fields as needed
            _speciesData[speciesId] = speciesData;
            
            std::cout << "Initialized species data for species " << speciesId 
                      << " with " << genomeIndices.size() << " genomes" << std::endl;
        }
    }

protected:
    // Fitness tracking components
    std::unique_ptr<Operator::PlotElitesParams> _eliteParams;
    std::unordered_map<uint32_t, DynamicSpeciesData> _speciesData;
    
    // Helper function to evaluate realistic fitness on evolved population
    void evaluateRealisticFitnessOnEvolvedPopulation() {
        auto& evolvedGenomes = _populationContainer->getCurrentGenomes(_foundationGeneration);
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(_foundationGeneration);
        
        // Clear any existing fitness results
        currentFitnessResults.clear();
        
        // Evaluate fitness for each evolved genome
        for (size_t i = 0; i < evolvedGenomes.size(); ++i) {
            const auto& genome = evolvedGenomes[i];
            
            // Use a simple fitness function based on genome complexity
            // This creates realistic fitness variation for testing elite selection
            double fitnessValue = calculateComplexityBasedFitness(genome);
            TestFitnessResult fitness(fitnessValue);
            
            currentFitnessResults.insert({fitness, i});
            
            std::cout << "Genome " << i << " fitness: " << fitnessValue << std::endl;
        }
        
        std::cout << "Evaluated fitness for " << evolvedGenomes.size() << " evolved genomes" << std::endl;
    }
    
    // Helper function to calculate fitness based on genome structural complexity
    double calculateComplexityBasedFitness(const Genome& genome) {
        const auto& nodeGenes = genome.get_nodeGenes();
        const auto& connectionGenes = genome.get_connectionGenes();
        
        // Simple fitness function: reward structural complexity
        double nodeScore = nodeGenes.size() * 10.0;
        double connectionScore = connectionGenes.size() * 5.0;
        
        // Add some variation based on node types
        double activationVariety = 0.0;
        std::set<ActivationType> uniqueActivations;
        for (const auto& node : nodeGenes) {
            uniqueActivations.insert(node.get_attributes().activationType);
        }
        activationVariety = uniqueActivations.size() * 3.0;
        
        // Add small random component for fitness variation within same structure
        std::hash<std::string> hasher;
        size_t genomeHash = hasher(std::to_string(nodeGenes.size()) + "_" + std::to_string(connectionGenes.size()));
        double randomComponent = (genomeHash % 100) / 100.0; // 0.0 to 0.99
        
        return nodeScore + connectionScore + activationVariety + randomComponent;
    }
    
    // Helper function to test elite selection on evolved species
    void testEliteSelectionOnEvolvedSpecies() {
        auto& evolvedGenomeData = _populationContainer->getCurrentGenomeData(_foundationGeneration);
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(_foundationGeneration);
        
        // Perform species grouping on the evolved population with fitness
        auto speciesGrouping = Population::speciesGrouping(
            currentFitnessResults, 
            evolvedGenomeData, 
            _speciesData, 
            *_globalIndexRegistry
        );
        
        std::cout << "\n=== SPECIES GROUPING RESULTS ===" << std::endl;
        for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
            std::cout << "Species " << speciesId << ": " << genomeIndices.size() << " genomes" << std::endl;
        }
        
        // Test elite selection algorithm
        std::cout << "\n=== BEFORE ELITE SELECTION ===" << std::endl;
        logEliteStatus();
        
        // Run PlotElites on the evolved population
        Operator::plotElites(speciesGrouping, *_eliteParams, *_globalIndexRegistry, _speciesData);
        
        std::cout << "\n=== AFTER ELITE SELECTION ===" << std::endl;
        logEliteStatus();
        
        // Validate elite selection results
        validateEliteSelectionResults(speciesGrouping, currentFitnessResults);
    }
    
    // Helper function to log current elite status for debugging
    void logEliteStatus() {
        auto& evolvedGenomes = _populationContainer->getCurrentGenomes(_foundationGeneration);
        auto& evolvedGenomeData = _populationContainer->getCurrentGenomeData(_foundationGeneration);
        
        for (size_t i = 0; i < evolvedGenomes.size(); ++i) {
            auto state = _globalIndexRegistry->getState(i);
            uint32_t speciesId = evolvedGenomeData[i].speciesId;
            bool isElite = (state == GenomeState::Elite);
            
            std::cout << "Genome " << i << ": species=" << speciesId 
                      << ", state=" << static_cast<int>(state) 
                      << ", elite=" << isElite << std::endl;
        }
    }
    
    // Helper function to validate that elite selection worked correctly
    void validateEliteSelectionResults(
        const std::unordered_map<uint32_t, std::vector<size_t>>& speciesGrouping,
        const std::multimap<TestFitnessResult, size_t>& fitnessResults) {
        
        // For each species, verify that at least one genome with the highest fitness is marked as elite
        for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
            // Find the highest fitness value in this species
            TestFitnessResult bestFitness(0.0);
            bool foundBestFitness = false;
            
            for (size_t genomeIndex : genomeIndices) {
                for (const auto& [fitness, globalIndex] : fitnessResults) {
                    if (globalIndex == genomeIndex) {
                        if (!foundBestFitness || fitness.isBetterThan(bestFitness)) {
                            bestFitness = fitness;
                            foundBestFitness = true;
                        }
                        break;
                    }
                }
            }
            
            // Find all genomes with the best fitness and verify at least one is elite
            if (foundBestFitness) {
                std::vector<size_t> bestGenomes;
                std::vector<size_t> bestEliteGenomes;
                
                for (size_t genomeIndex : genomeIndices) {
                    for (const auto& [fitness, globalIndex] : fitnessResults) {
                        if (globalIndex == genomeIndex) {
                            if (fitness.isEqualTo(bestFitness)) {
                                bestGenomes.push_back(genomeIndex);
                                auto state = _globalIndexRegistry->getState(genomeIndex);
                                if (state == GenomeState::Elite) {
                                    bestEliteGenomes.push_back(genomeIndex);
                                }
                            }
                            break;
                        }
                    }
                }
                
                // Verify at least one genome with best fitness is marked as elite
                assert(!bestEliteGenomes.empty() && 
                       "At least one genome with best fitness in species should be marked as Elite");
                
                std::cout << "Species " << speciesId << ": " << bestEliteGenomes.size() 
                          << " elite(s) among " << bestGenomes.size() 
                          << " genome(s) with best fitness " << bestFitness.getFitness() << std::endl;
            }
        }
    }
};

// ============================================================================
// ELITE PRESERVATION TEST FIXTURE - Level 3 Integration Tests  
// ============================================================================

// Level 3 fixture that inherits fitness tracking from Level 2 and adds multi-generation
// elite preservation testing. This tests the exact mechanisms that cause the species
// retention bug by validating elite preservation across multiple generations.
class ElitePreservationTest : public FitnessTrackingTest {
protected:
    void SetUp() override {
        // Inherit all setup from Level 2 (foundation + fitness tracking + elite selection)
        FitnessTrackingTest::SetUp();
        
        // Perform initial elite selection to establish baseline
        establishEliteBaseline();
        
        std::cout << "Level 3 Elite Preservation Test ready with elite baseline established" << std::endl;
    }

private:
    void establishEliteBaseline() {
        // Run initial fitness evaluation and elite selection (from Level 2)
        evaluateRealisticFitnessOnEvolvedPopulation();
        testEliteSelectionOnEvolvedSpecies();
        
        // Capture elite baseline state for multi-generation tracking
        captureEliteBaseline();
        
        std::cout << "Elite baseline established: " << _baselineElites.size() 
                  << " elites across " << _baselineElitesBySpecies.size() << " species" << std::endl;
    }
    
    void captureEliteBaseline() {
        auto& evolvedGenomes = _populationContainer->getCurrentGenomes(_foundationGeneration);
        auto& evolvedGenomeData = _populationContainer->getCurrentGenomeData(_foundationGeneration);
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(_foundationGeneration);
        
        // Clear baseline data
        _baselineElites.clear();
        _baselineEliteProperties.clear();
        _baselineEliteFitness.clear();
        _baselineEliteSpecies.clear();
        _baselineElitesBySpecies.clear();
        
        // Identify all current elites
        for (size_t i = 0; i < evolvedGenomes.size(); ++i) {
            auto state = _globalIndexRegistry->getState(i);
            if (state == GenomeState::Elite) {
                _baselineElites.push_back(i);
                
                // Capture elite properties for byte-for-byte comparison
                GenomeProperties props = extractGenomeProperties(evolvedGenomes[i]);
                _baselineEliteProperties.push_back(props);
                
                // Capture elite fitness
                TestFitnessResult eliteFitness(0.0);
                for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                    if (globalIndex == i) {
                        eliteFitness = fitness;
                        break;
                    }
                }
                _baselineEliteFitness.push_back(eliteFitness);
                
                // Capture elite species
                uint32_t speciesId = evolvedGenomeData[i].speciesId;
                _baselineEliteSpecies.push_back(speciesId);
                _baselineElitesBySpecies[speciesId].push_back(i);
                
                std::cout << "Elite baseline: genome " << i << " in species " << speciesId 
                          << " with fitness " << eliteFitness.getFitness() << std::endl;
            }
        }
    }

protected:
    // Elite baseline data for multi-generation tracking
    std::vector<size_t> _baselineElites;                                    // Elite genome indices
    std::vector<GenomeProperties> _baselineEliteProperties;                 // Elite genome properties for comparison
    std::vector<TestFitnessResult> _baselineEliteFitness;                   // Elite fitness values
    std::vector<uint32_t> _baselineEliteSpecies;                           // Elite species assignments
    std::unordered_map<uint32_t, std::vector<size_t>> _baselineElitesBySpecies; // Elites grouped by species
    
    // Helper function to run a single generation of evolution with elite protection
    void runSingleGenerationEvolution(uint32_t targetGeneration) {
        // This simulates the 1:1 evolution phase from EvolutionPrototype
        // with elite protection logic active
        
        auto& lastGenomes = _populationContainer->getLastGenomes(targetGeneration);
        auto& lastGenomeData = _populationContainer->getLastGenomeData(targetGeneration);
        auto& currentGenomes = _populationContainer->getCurrentGenomes(targetGeneration);
        auto& currentGenomeData = _populationContainer->getCurrentGenomeData(targetGeneration);
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(targetGeneration);
        
        // Clear current generation fitness results
        currentFitnessResults.clear();
        
        std::cout << "\n=== RUNNING GENERATION " << targetGeneration << " EVOLUTION ===" << std::endl;
        std::cout << "Processing " << lastGenomes.size() << " genomes from previous generation" << std::endl;
        
        // Phase 1: 1:1 Evolution with Elite Protection (similar to EvolutionPrototype)
        for (size_t i = 0; i < lastGenomes.size(); ++i) {
            const Genome& parentGenome = lastGenomes[i];
            const Population::DynamicGenomeData& parentData = lastGenomeData[i];
            
            // Elite protection logic (exact copy from EvolutionPrototype.hpp:441-448)
            if (_globalIndexRegistry->getState(i) == GenomeState::Elite) {
                // Elite protection: copy as-is without any mutation
                currentGenomes[i] = parentGenome;
                currentGenomeData[i] = parentData;
                
                // Verify elite copy is identical (assert from EvolutionPrototype)
                bool isIdentical = Operator::genomeEquals(parentGenome, currentGenomes[i]);
                EXPECT_TRUE(isIdentical) << "Elite genome " << i << " copy should be identical to parent";
                
                std::cout << "Elite " << i << " protected and copied identically" << std::endl;
            } else {
                // Non-elite genomes: apply basic mutation (simplified for testing)
                currentGenomes[i] = applySimpleMutation(parentGenome);
                currentGenomeData[i] = parentData;
                
                std::cout << "Non-elite " << i << " mutated" << std::endl;
            }
        }
        
        // Phase 2: Fitness Evaluation with Elite Fitness Preservation (EvolutionPrototype.hpp:544-565)
        for (size_t i = 0; i < currentGenomes.size(); ++i) {
            const Genome& genome = currentGenomes[i];
            
            if (_globalIndexRegistry->getState(i) == GenomeState::Elite) {
                // Elite fitness preservation: copy fitness from previous generation
                const auto& lastFitnessResults = _populationContainer->getLastFitnessResults(targetGeneration);
                TestFitnessResult eliteFitness(0.0);
                bool fitnessFound = false;
                
                for (const auto& [fitness, globalIndex] : lastFitnessResults) {
                    if (globalIndex == i) {
                        eliteFitness = fitness;
                        fitnessFound = true;
                        break;
                    }
                }
                
                EXPECT_TRUE(fitnessFound) << "Elite genome " << i << " should have fitness from previous generation";
                
                // Preserve elite fitness (no re-evaluation)
                currentFitnessResults.insert({eliteFitness, i});
                
                std::cout << "Elite " << i << " fitness preserved: " << eliteFitness.getFitness() << std::endl;
            } else {
                // Non-elite fitness evaluation
                double fitnessValue = calculateComplexityBasedFitness(genome);
                TestFitnessResult fitness(fitnessValue);
                currentFitnessResults.insert({fitness, i});
                
                // Update species assignment for non-elites
                uint32_t newSpeciesId = Operator::compatibilityDistance(genome, _historyTracker, *_compatibilityParams);
                currentGenomeData[i].speciesId = newSpeciesId;
            }
        }
        
        std::cout << "Generation " << targetGeneration << " evolution completed" << std::endl;
    }
    
    // Helper function to apply simple mutation for non-elite genomes
    Genome applySimpleMutation(const Genome& parent) {
        // Apply a basic weight mutation to create variation
        Operator::WeightMutationParams weightParams(0.8, 0.1, 0.1, 2.0, Operator::WeightMutationParams::MutationType::MIXED);
        Genome offspring = Operator::weightMutation(parent, weightParams);
        Operator::phenotypeUpdateWeight(offspring);
        return offspring;
    }
    
    // Helper function to validate elite identity preservation
    void validateEliteIdentityPreservation(uint32_t generation) {
        auto& currentGenomes = _populationContainer->getCurrentGenomes(generation);
        
        std::cout << "\n=== VALIDATING ELITE IDENTITY PRESERVATION ===" << std::endl;
        
        for (size_t i = 0; i < _baselineElites.size(); ++i) {
            size_t eliteIndex = _baselineElites[i];
            
            // Extract current properties
            GenomeProperties currentProps = extractGenomeProperties(currentGenomes[eliteIndex]);
            
            // Compare with baseline (must be identical)
            EXPECT_TRUE(currentProps.isIdenticalTo(_baselineEliteProperties[i]))
                << "Elite genome " << eliteIndex << " properties changed in generation " << generation;
                
            std::cout << "Elite " << eliteIndex << " identity preserved: " 
                      << (currentProps.isIdenticalTo(_baselineEliteProperties[i]) ? "PASS" : "FAIL") << std::endl;
        }
    }
    
    // Helper function to validate elite fitness preservation
    void validateEliteFitnessPreservation(uint32_t generation) {
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(generation);
        
        std::cout << "\n=== VALIDATING ELITE FITNESS PRESERVATION ===" << std::endl;
        
        for (size_t i = 0; i < _baselineElites.size(); ++i) {
            size_t eliteIndex = _baselineElites[i];
            
            // Find current fitness
            TestFitnessResult currentFitness(0.0);
            bool found = false;
            for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                if (globalIndex == eliteIndex) {
                    currentFitness = fitness;
                    found = true;
                    break;
                }
            }
            
            EXPECT_TRUE(found) << "Elite genome " << eliteIndex << " should have fitness in generation " << generation;
            
            if (found) {
                // Compare with baseline fitness (must be identical)
                EXPECT_TRUE(currentFitness.isEqualTo(_baselineEliteFitness[i]))
                    << "Elite genome " << eliteIndex << " fitness changed from " 
                    << _baselineEliteFitness[i].getFitness() << " to " << currentFitness.getFitness()
                    << " in generation " << generation;
                    
                std::cout << "Elite " << eliteIndex << " fitness preserved: " 
                          << currentFitness.getFitness() << " == " << _baselineEliteFitness[i].getFitness()
                          << " (" << (currentFitness.isEqualTo(_baselineEliteFitness[i]) ? "PASS" : "FAIL") << ")" << std::endl;
            }
        }
    }
    
    // Helper function to validate species fitness progression (the critical assertion!)
    void validateSpeciesFitnessProgression(uint32_t generation) {
        auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(generation);
        auto& currentGenomeData = _populationContainer->getCurrentGenomeData(generation);
        
        std::cout << "\n=== VALIDATING SPECIES FITNESS PROGRESSION (CRITICAL ASSERTION) ===" << std::endl;
        
        // Test the exact logic from EvolutionPrototype.hpp:706-733
        for (const auto& [speciesId, eliteIndices] : _baselineElitesBySpecies) {
            // Find current best fitness for this species
            TestFitnessResult currentBest(0.0);
            bool foundGenome = false;
            
            for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                if (currentGenomeData[globalIndex].speciesId == speciesId) {
                    auto state = _globalIndexRegistry->getState(globalIndex);
                    if ((state == GenomeState::Active || state == GenomeState::Elite)) {
                        if (!foundGenome || fitness.isBetterThan(currentBest)) {
                            currentBest = fitness;
                            foundGenome = true;
                        }
                    }
                }
            }
            
            if (foundGenome) {
                // Find previous best fitness for this species (from elite baseline)
                TestFitnessResult previousBest(0.0);
                for (size_t eliteIdx : eliteIndices) {
                    auto eliteIt = std::find(_baselineElites.begin(), _baselineElites.end(), eliteIdx);
                    if (eliteIt != _baselineElites.end()) {
                        size_t baselineIndex = std::distance(_baselineElites.begin(), eliteIt);
                        if (previousBest.getFitness() == 0.0 || _baselineEliteFitness[baselineIndex].isBetterThan(previousBest)) {
                            previousBest = _baselineEliteFitness[baselineIndex];
                        }
                    }
                }
                
                // This is the EXACT assertion that fails in EvolutionPrototype!
                bool progressionValid = currentBest.isBetterThan(previousBest) || currentBest.isEqualTo(previousBest);
                EXPECT_TRUE(progressionValid)
                    << "SPECIES FITNESS REGRESSION in generation " << generation 
                    << " for species " << speciesId
                    << ": current=" << currentBest.getFitness() 
                    << " vs previous=" << previousBest.getFitness()
                    << " (assertion from EvolutionPrototype.hpp:731-732)";
                    
                std::cout << "Species " << speciesId << " fitness progression: " 
                          << previousBest.getFitness() << " -> " << currentBest.getFitness()
                          << " (" << (progressionValid ? "PASS" : "FAIL") << ")" << std::endl;
            }
        }
    }
    
    // Helper function to validate elite status retention
    void validateEliteStatusRetention(uint32_t generation) {
        std::cout << "\n=== VALIDATING ELITE STATUS RETENTION ===" << std::endl;
        
        for (size_t eliteIndex : _baselineElites) {
            auto currentState = _globalIndexRegistry->getState(eliteIndex);
            bool isStillElite = (currentState == GenomeState::Elite);
            
            EXPECT_TRUE(isStillElite) 
                << "Elite genome " << eliteIndex << " lost Elite status in generation " << generation
                << " (current state: " << static_cast<int>(currentState) << ")";
                
            std::cout << "Elite " << eliteIndex << " status retention: " 
                      << (isStillElite ? "PASS" : "FAIL") << std::endl;
        }
    }
};

// Level 1: Core Genome Propagation Integration Test
// Tests coordination between PopulationContainer, GlobalIndexRegistry, 
// HistoryTracker, and CompatibilityDistance through proper evolutionary process
TEST_F(CoreGenomePropagationTest, EvolutionaryProcess_GenomeConsistency) {
    // ========================================================================
    // SIMPLIFIED TEST - Debug the genome copying issue
    // ========================================================================
    
    const size_t INITIAL_POPULATION = 7;  // Start small to debug
    const size_t GENERATIONS = 8;        // Only test a few generations
    const size_t NEW_PER_GENERATION = 5;  // Create 3 new genomes through mutation each generation
    
    // Create Init operator parameters for consistent genome creation
    std::vector<NodeGeneAttributes> inputAttributes = {{ActivationType::NONE}};  // 1 input only
    std::vector<NodeGeneAttributes> outputAttributes = {{ActivationType::SIGMOID}};  // 1 output
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttributes;
    biasAttributes[0] = {1.0f, true};  // Bias to output connection
    
    Operator::InitParams initParams(inputAttributes, outputAttributes, biasAttributes, 
                                  Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    // Store baseline properties for all genomes across generations
    std::vector<std::vector<GenomeProperties>> baselineProperties(GENERATIONS);
    
    // Store species assignments for validation
    std::vector<std::vector<uint32_t>> baselineSpeciesIDs(GENERATIONS);
    
    // Mutation probability parameters (similar to EvolutionPrototype)
    MutationProbabilities mutationProbs;
    mutationProbs.weightMutation = 0.6;      // High probability for weight mutations
    mutationProbs.nodeMutation = 0.15;       // Low probability for structural changes  
    mutationProbs.connectionMutation = 0.2;  // Low probability for structural changes
    mutationProbs.connectionReactivation = 0.05; // Very low probability for reactivation
    
    // Mutation operator parameters
    Operator::WeightMutationParams weightMutationParams(0.8, 0.1, 0.1, 2.0, Operator::WeightMutationParams::MutationType::MIXED);
    Operator::NodeMutationParams nodeMutationParams({ActivationType::TANH});  // New hidden nodes use TANH
    Operator::ConnectionMutationParams connectionMutationParams(1.0, 2.0, Operator::ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    Operator::ConnectionReactivationParams connectionReactivationParams(Operator::ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    // Random number generator for mutations
    std::random_device rd;
    std::mt19937 gen_rand(rd());
    
    // ========================================================================
    // GENERATION 0 - Create Initial Population Using Init Operator
    // ========================================================================
    
    baselineProperties[0].reserve(INITIAL_POPULATION);
    baselineSpeciesIDs[0].reserve(INITIAL_POPULATION);
    
    // Create initial population using Init operator
    for (size_t i = 0; i < INITIAL_POPULATION; ++i) {
        // Create genome using Init operator
        Genome genome = Operator::init(_historyTracker, initParams);
        
        // Test: Extract properties immediately after creation
        GenomeProperties props = extractGenomeProperties(genome);
        
        // Calculate species assignment using compatibility distance
        uint32_t speciesId = Operator::compatibilityDistance(genome, _historyTracker, *_compatibilityParams);
        
        // Create genome data with proper species assignment
        Population::DynamicGenomeData genomeData;
        genomeData.speciesId = speciesId;
        genomeData.pendingEliminationCounter = 0;
        genomeData.repairAttempts = 0;
        genomeData.isUnderRepair = false;
        genomeData.isMarkedForElimination = false;
        genomeData.genomeIndex = UINT32_MAX;  // Will be set by container
        genomeData.parentAIndex = UINT32_MAX;  // No parents for initial genomes
        genomeData.parentBIndex = UINT32_MAX;
        
        // Store baseline properties and species ID
        baselineProperties[0].push_back(props);
        baselineSpeciesIDs[0].push_back(speciesId);
        
        std::cout << "Gen 0 genome " << i << " assigned to species " << speciesId << std::endl;
        
        // Add to population
        uint32_t index = _populationContainer->push_back(0, std::move(genome), genomeData);
        EXPECT_EQ(index, i) << "Initial genome index mismatch at position " << i;
    }
    
    // Validate generation 0 setup
    std::cout << "Generation 0 created " << _populationContainer->getGenerationSize(0) << " genomes" << std::endl;
    EXPECT_EQ(_populationContainer->getGenerationSize(0), INITIAL_POPULATION);
    validateGlobalIndexConsistency(*_populationContainer, *_globalIndexRegistry, 0);
    validateHistoryTrackerConsistency(_historyTracker, 0);
    validateSpeciesConsistency(0, baselineSpeciesIDs, 0); // No copied genomes in generation 0
    
    // ========================================================================
    // SIMPLIFIED EVOLUTIONARY LOOP - Test basic copying
    // ========================================================================
    
    for (uint32_t gen = 1; gen < GENERATIONS; ++gen) {
        // Get the actual size of the previous generation
        size_t prevGenSize = _populationContainer->getGenerationSize(gen - 1);
        std::cout << "Gen " << gen << " - Previous generation actual size: " << prevGenSize << std::endl;
        
        if (prevGenSize == 0) {
            std::cout << "Gen " << gen << " - No genomes in previous generation to copy, skipping" << std::endl;
            continue;
        }
        
        // Reserve space for copied genomes + new mutated genomes
        baselineProperties[gen].reserve(prevGenSize + NEW_PER_GENERATION);
        baselineSpeciesIDs[gen].reserve(prevGenSize + NEW_PER_GENERATION);
        
        // Step 1: DIRECT ASSIGNMENT - Copy ALL genomes from last generation to current generation
        auto& lastGenomes = _populationContainer->getLastGenomes(gen);
        auto& lastGenomeData = _populationContainer->getLastGenomeData(gen);
        auto& currentGenomes = _populationContainer->getCurrentGenomes(gen);
        auto& currentGenomeData = _populationContainer->getCurrentGenomeData(gen);
        
        std::cout << "Gen " << gen << " - Direct copying " << prevGenSize << " genomes from last generation" << std::endl;
        
        // Direct assignment of all genomes and their data
        currentGenomes = lastGenomes;
        currentGenomeData = lastGenomeData;
        
        // Store baselines for all copied genomes
        for (size_t i = 0; i < prevGenSize; ++i) {
            // Copy baseline from previous generation
            if (i < baselineProperties[gen - 1].size()) {
                baselineProperties[gen].push_back(baselineProperties[gen - 1][i]);
                baselineSpeciesIDs[gen].push_back(baselineSpeciesIDs[gen - 1][i]);
            }
            
            std::cout << "Gen " << gen << " copied genome " << i 
                      << " with " << currentGenomes[i].get_nodeGenes().size() << " nodes, "
                      << currentGenomes[i].get_connectionGenes().size() << " connections" << std::endl;
        }
        
        // ========================================================================
        // Step 2: CREATE NEW GENOMES THROUGH MUTATIONS (using last NEW_PER_GENERATION genomes as parents)
        // ========================================================================
        
        std::cout << "Gen " << gen << " - Creating " << NEW_PER_GENERATION << " new genomes through mutations" << std::endl;
        
        // Calculate the range of last NEW_PER_GENERATION genomes to use as parents
        size_t startParentIdx = (prevGenSize >= NEW_PER_GENERATION) ? (prevGenSize - NEW_PER_GENERATION) : 0;
        std::cout << "Gen " << gen << " - Using parents from index " << startParentIdx << " to " << (prevGenSize - 1) << std::endl;
        
        for (size_t newIdx = 0; newIdx < NEW_PER_GENERATION; ++newIdx) {
            // Select parent from the last NEW_PER_GENERATION genomes
            size_t parentIdx = startParentIdx + newIdx;
            if (parentIdx >= prevGenSize) {
                parentIdx = prevGenSize - 1; // Fallback to last genome if needed
            }
            
            const auto& parentGenome = lastGenomes[parentIdx];
            std::cout << "Gen " << gen << " - Creating offspring " << newIdx << " from parent " << parentIdx << std::endl;
            
            // Apply mutations with connection validation
            Genome offspring = applyMutationWithValidation(
                parentGenome, 
                _historyTracker, 
                mutationProbs,
                weightMutationParams,
                nodeMutationParams,
                connectionMutationParams,
                connectionReactivationParams,
                gen_rand);
            
            // Calculate species assignment for the mutated genome
            uint32_t newSpeciesId = Operator::compatibilityDistance(offspring, _historyTracker, *_compatibilityParams);
            uint32_t parentSpeciesId = baselineSpeciesIDs[gen][parentIdx]; // Parent is now in current generation
            
            std::cout << "Gen " << gen << " - Offspring " << newIdx 
                      << " with " << offspring.get_nodeGenes().size() << " nodes, "
                      << offspring.get_connectionGenes().size() << " connections"
                      << " (parent species " << parentSpeciesId << " -> offspring species " << newSpeciesId << ")" << std::endl;
            
            // Create genome data for new offspring
            Population::DynamicGenomeData newGenomeData;
            newGenomeData.speciesId = newSpeciesId;
            newGenomeData.pendingEliminationCounter = 0;
            newGenomeData.repairAttempts = 0;
            newGenomeData.isUnderRepair = false;
            newGenomeData.isMarkedForElimination = false;
            newGenomeData.genomeIndex = UINT32_MAX;  // Will be set by container
            newGenomeData.parentAIndex = parentIdx;     // Track lineage
            newGenomeData.parentBIndex = UINT32_MAX;    // Single parent (mutation, not crossover)
            
            // Store baseline properties for new genome
            GenomeProperties newProps = extractGenomeProperties(offspring);
            baselineProperties[gen].push_back(newProps);
            baselineSpeciesIDs[gen].push_back(newSpeciesId);
            
            // Add new genome to population
            _populationContainer->push_back(gen, std::move(offspring), newGenomeData);
        }
        
        // === CRITICAL VALIDATIONS ===
        
        // 1. Global indices remain valid and accessible
        validateGlobalIndexConsistency(*_populationContainer, *_globalIndexRegistry, gen);
        
        // 2. History tracker maintains innovation number consistency
        validateHistoryTrackerConsistency(_historyTracker, gen);
        
        // 3. Species assignment consistency (only for copied genomes)
        validateSpeciesConsistency(gen, baselineSpeciesIDs, prevGenSize);
        
        // 4. Genome content consistency (only for copied genomes)
        validateAllGenomeConsistency(gen, baselineProperties, prevGenSize);
        
        // 5. Population size validation (accounts for PopulationContainer synchronization)
        // Each generation adds prevGenSize (copied) + NEW_PER_GENERATION (new) genomes
        size_t expectedSize = prevGenSize + NEW_PER_GENERATION;
        EXPECT_EQ(_populationContainer->getGenerationSize(gen), expectedSize) 
            << "Population size incorrect in generation " << gen;
    }
    
    // ========================================================================
    // FINAL VALIDATION - System State After Evolutionary Process
    // ========================================================================
    
    // Final population size calculation
    // Gen 0: 3, Gen 1: 3+3=6, Gen 2: 6+3=9, etc.
    // Each generation copies from previous + adds NEW_PER_GENERATION new genomes
    size_t expectedFinalSize = baselineProperties[GENERATIONS - 1].size();
    EXPECT_EQ(_populationContainer->getGenerationSize(GENERATIONS - 1), expectedFinalSize);
    
    // Registry should accommodate all genomes
    EXPECT_GE(_globalIndexRegistry->getMaxIndex(), expectedFinalSize);
    
    // All generations should maintain consistent sizing due to synchronization
    EXPECT_TRUE(_populationContainer->validateConsistency());
    
    // Final comprehensive validation of all systems
    size_t finalCopiedGenomes = baselineProperties[GENERATIONS - 2].size(); // Previous generation size
    validateSpeciesConsistency(GENERATIONS - 1, baselineSpeciesIDs, finalCopiedGenomes);
    validateAllGenomeConsistency(GENERATIONS - 1, baselineProperties, finalCopiedGenomes);
    
    // Summary validation: Verify species diversity and consistency across all generations
    std::cout << "\n=== SPECIES ASSIGNMENT SUMMARY ===" << std::endl;
    for (uint32_t gen = 0; gen < GENERATIONS; ++gen) {
        std::cout << "Generation " << gen << " species assignments: ";
        for (size_t i = 0; i < baselineSpeciesIDs[gen].size(); ++i) {
            std::cout << baselineSpeciesIDs[gen][i] << " ";
        }
        std::cout << std::endl;
    }
    
    // Innovation tracking validation - verify structural diversity has increased
    std::cout << "\n=== INNOVATION TRACKING VALIDATION ===" << std::endl;
    
    // Count unique node and connection history IDs across all final genomes
    std::set<uint32_t> allNodeHistoryIDs;
    std::set<uint32_t> allConnectionHistoryIDs;
    
    auto& finalGenomes = _populationContainer->getCurrentGenomes(GENERATIONS - 1);
    for (const auto& genome : finalGenomes) {
        for (const auto& node : genome.get_nodeGenes()) {
            allNodeHistoryIDs.insert(node.get_historyID());
        }
        for (const auto& conn : genome.get_connectionGenes()) {
            allConnectionHistoryIDs.insert(conn.get_historyID());
        }
    }
    
    std::cout << "Final population structural diversity:" << std::endl;
    std::cout << "  Unique node innovation IDs: " << allNodeHistoryIDs.size() << std::endl;
    std::cout << "  Unique connection innovation IDs: " << allConnectionHistoryIDs.size() << std::endl;
    
    // Verify that structural evolution occurred (more innovations than initial)
    // Initial genomes have: 1 input + 1 bias + 1 output = 3 nodes, 2 connections (input->output, bias->output)
    EXPECT_GE(allNodeHistoryIDs.size(), 3u) << "Should have at least initial node innovations";
    EXPECT_GE(allConnectionHistoryIDs.size(), 2u) << "Should have at least initial connection innovations";
    
    // If mutations were applied successfully, we should see structural growth
    if (allNodeHistoryIDs.size() > 3) {
        std::cout << "  SUCCESS: Structural evolution detected - new nodes were added!" << std::endl;
    }
    if (allConnectionHistoryIDs.size() > 2) {
        std::cout << "  SUCCESS: Structural evolution detected - new connections were added!" << std::endl;
    }
    
    // SUCCESS: If we reach here, genome copying and species management work correctly
}

// ============================================================================
// LEVEL 2: ELITE SELECTION AND FITNESS CORRELATION INTEGRATION TESTS
// ============================================================================

// Test 2: Elite Selection Accuracy - Tests PlotElites on realistically evolved population
// This test inherits the evolved foundation from Level 1 and validates elite selection
TEST_F(FitnessTrackingTest, EliteSelectionAccuracy) {
    // ========================================================================
    // PHASE 1: FITNESS EVALUATION - On Already Evolved Population
    // ========================================================================
    
    std::cout << "\n=== PHASE 1: FITNESS EVALUATION ===" << std::endl;
    
    // Evaluate fitness on the realistic evolved population from Level 1
    evaluateRealisticFitnessOnEvolvedPopulation();
    
    // Access the evolved state (provided by EvolutionaryFoundationFixture)
    auto& evolvedGenomes = _populationContainer->getCurrentGenomes(_foundationGeneration);
    auto& evolvedGenomeData = _populationContainer->getCurrentGenomeData(_foundationGeneration);
    auto& currentFitnessResults = _populationContainer->getCurrentFitnessResults(_foundationGeneration);
    
    // Validate that we have meaningful evolved state
    EXPECT_GT(_evolvedPopulationSize, 0u) << "Should have evolved population";
    EXPECT_GT(_evolvedSpeciesMap.size(), 0u) << "Should have species diversity";
    EXPECT_EQ(evolvedGenomes.size(), _evolvedPopulationSize) << "Genome count should match evolved population size";
    EXPECT_EQ(currentFitnessResults.size(), _evolvedPopulationSize) << "Should have fitness for every genome";
    
    std::cout << "Evolved population validated: " << _evolvedPopulationSize 
              << " genomes across " << _evolvedSpeciesMap.size() << " species" << std::endl;
    
    // ========================================================================
    // PHASE 2: SPECIES GROUPING VALIDATION - Verify Fitness Ordering
    // ========================================================================
    
    std::cout << "\n=== PHASE 2: SPECIES GROUPING VALIDATION ===" << std::endl;
    
    // Perform species grouping on the evolved population with fitness
    auto speciesGrouping = Population::speciesGrouping(
        currentFitnessResults, 
        evolvedGenomeData, 
        _speciesData, 
        *_globalIndexRegistry
    );
    
    // Validate species grouping structure
    EXPECT_GT(speciesGrouping.size(), 0u) << "Should have at least one species group";
    EXPECT_LE(speciesGrouping.size(), _evolvedSpeciesMap.size()) 
        << "Species grouping should not exceed evolved species count";
    
    // Validate fitness ordering within each species (worst to best)
    for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
        EXPECT_GT(genomeIndices.size(), 0u) << "Species " << speciesId << " should have genomes";
        
        std::cout << "Species " << speciesId << " fitness ordering: ";
        std::vector<double> speciesFitnessValues;
        
        for (size_t genomeIndex : genomeIndices) {
            // Find fitness for this genome
            for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                if (globalIndex == genomeIndex) {
                    double fitnessValue = fitness.getFitness();
                    speciesFitnessValues.push_back(fitnessValue);
                    std::cout << fitnessValue << " ";
                    break;
                }
            }
        }
        std::cout << std::endl;
        
        // Verify ascending fitness order (worst to best)
        for (size_t i = 1; i < speciesFitnessValues.size(); ++i) {
            EXPECT_GE(speciesFitnessValues[i], speciesFitnessValues[i-1])
                << "Fitness should be ordered worst to best in species " << speciesId;
        }
    }
    
    // ========================================================================
    // PHASE 3: ELITE SELECTION TESTING - Test PlotElites Algorithm
    // ========================================================================
    
    std::cout << "\n=== PHASE 3: ELITE SELECTION TESTING ===" << std::endl;
    
    // Log status before elite selection
    std::cout << "Before elite selection:" << std::endl;
    logEliteStatus();
    
    // Calculate expected elite count per species
    std::unordered_map<uint32_t, size_t> expectedElitesBySpecies;
    for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
        size_t speciesSize = genomeIndices.size();
        size_t expectedElites = static_cast<size_t>(std::ceil(speciesSize * 0.20)); // 20%
        expectedElites = std::max(expectedElites, static_cast<size_t>(1)); // min 1
        expectedElites = std::min(expectedElites, static_cast<size_t>(3)); // max 3
        expectedElites = std::min(expectedElites, speciesSize); // can't exceed species size
        
        expectedElitesBySpecies[speciesId] = expectedElites;
        std::cout << "Species " << speciesId << ": " << speciesSize << " genomes, expecting " 
                  << expectedElites << " elites" << std::endl;
    }
    
    // Run PlotElites algorithm on the evolved population
    Operator::plotElites(speciesGrouping, *_eliteParams, *_globalIndexRegistry, _speciesData);
    
    // Log status after elite selection
    std::cout << "After elite selection:" << std::endl;
    logEliteStatus();
    
    // ========================================================================
    // PHASE 4: ELITE SELECTION VALIDATION - Verify Correctness
    // ========================================================================
    
    std::cout << "\n=== PHASE 4: ELITE SELECTION VALIDATION ===" << std::endl;
    
    // Count actual elites per species
    std::unordered_map<uint32_t, size_t> actualElitesBySpecies;
    std::vector<size_t> allEliteIndices;
    
    for (size_t i = 0; i < evolvedGenomes.size(); ++i) {
        auto state = _globalIndexRegistry->getState(i);
        if (state == GenomeState::Elite) {
            uint32_t speciesId = evolvedGenomeData[i].speciesId;
            actualElitesBySpecies[speciesId]++;
            allEliteIndices.push_back(i);
        }
    }
    
    // Validate elite count per species
    for (const auto& [speciesId, expectedCount] : expectedElitesBySpecies) {
        size_t actualCount = actualElitesBySpecies[speciesId];
        EXPECT_EQ(actualCount, expectedCount) 
            << "Species " << speciesId << " should have " << expectedCount 
            << " elites but got " << actualCount;
    }
    
    // Validate that elites have the highest fitness in their species
    for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
        // Find all elites in this species
        std::vector<size_t> speciesElites;
        for (size_t genomeIndex : genomeIndices) {
            auto state = _globalIndexRegistry->getState(genomeIndex);
            if (state == GenomeState::Elite) {
                speciesElites.push_back(genomeIndex);
            }
        }
        
        // For each elite, verify it has higher fitness than all non-elites in the species
        for (size_t eliteIndex : speciesElites) {
            TestFitnessResult eliteFitness(0.0);
            bool foundEliteFitness = false;
            
            // Find elite's fitness
            for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                if (globalIndex == eliteIndex) {
                    eliteFitness = fitness;
                    foundEliteFitness = true;
                    break;
                }
            }
            
            ASSERT_TRUE(foundEliteFitness) << "Elite genome " << eliteIndex << " should have fitness";
            
            // Check that elite has higher fitness than all non-elites in species
            for (size_t genomeIndex : genomeIndices) {
                auto state = _globalIndexRegistry->getState(genomeIndex);
                if (state != GenomeState::Elite) {
                    // Find non-elite's fitness
                    for (const auto& [fitness, globalIndex] : currentFitnessResults) {
                        if (globalIndex == genomeIndex) {
                            EXPECT_TRUE(eliteFitness.isBetterThan(fitness) || eliteFitness.isEqualTo(fitness))
                                << "Elite genome " << eliteIndex << " (fitness " << eliteFitness.getFitness()
                                << ") should have >= fitness than non-elite " << genomeIndex 
                                << " (fitness " << fitness.getFitness() << ") in species " << speciesId;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // ========================================================================
    // PHASE 5: SUMMARY AND SUCCESS VALIDATION
    // ========================================================================
    
    std::cout << "\n=== PHASE 5: ELITE SELECTION SUMMARY ===" << std::endl;
    
    size_t totalElites = allEliteIndices.size();
    size_t totalGenomes = evolvedGenomes.size();
    double elitePercentage = (totalElites * 100.0) / totalGenomes;
    
    std::cout << "Elite selection results:" << std::endl;
    std::cout << "  Total genomes: " << totalGenomes << std::endl;
    std::cout << "  Total elites: " << totalElites << std::endl;
    std::cout << "  Elite percentage: " << elitePercentage << "%" << std::endl;
    std::cout << "  Species with elites: " << actualElitesBySpecies.size() << "/" << speciesGrouping.size() << std::endl;
    
    // Validate overall elite selection success
    EXPECT_GT(totalElites, 0u) << "Should have selected at least one elite";
    EXPECT_LE(totalElites, totalGenomes) << "Cannot have more elites than genomes";
    
    // Validate that every species with genomes has at least one elite
    for (const auto& [speciesId, genomeIndices] : speciesGrouping) {
        EXPECT_GT(actualElitesBySpecies[speciesId], 0u) 
            << "Species " << speciesId << " should have at least one elite";
    }
    
    std::cout << "SUCCESS: Elite selection correctly identified highest-fitness genomes within species boundaries!" << std::endl;
    
    // This test validates that PlotElites correctly:
    // ✅ Identifies highest-fitness genomes in each species
    // ✅ Respects species boundaries during selection
    // ✅ Calculates elite percentages correctly with min/max constraints
    // ✅ Updates GlobalIndexRegistry with correct elite status
    // ✅ Works on realistically evolved populations with natural species diversity
}

// ============================================================================
// LEVEL 3: ELITE PRESERVATION ACROSS GENERATIONS INTEGRATION TESTS
// ============================================================================

// Test 3: Multi-Generation Elite Preservation - Tests the exact mechanisms causing species retention bug
// This test validates elite preservation across multiple generations and tests the failing assertion
TEST_F(ElitePreservationTest, MultiGenerationElitePreservation) {
    // ========================================================================
    // PHASE 1: ELITE BASELINE VALIDATION - Verify Setup from Level 2
    // ========================================================================
    
    std::cout << "\n=== PHASE 1: ELITE BASELINE VALIDATION ===" << std::endl;
    
    // Validate that we have a meaningful elite baseline from Level 2
    EXPECT_GT(_baselineElites.size(), 0u) << "Should have elites from Level 2 setup";
    EXPECT_GT(_baselineElitesBySpecies.size(), 0u) << "Should have elites across multiple species";
    EXPECT_EQ(_baselineElites.size(), _baselineEliteProperties.size()) << "Elite data consistency check";
    EXPECT_EQ(_baselineElites.size(), _baselineEliteFitness.size()) << "Elite fitness data consistency";
    
    std::cout << "Elite baseline validated: " << _baselineElites.size() 
              << " elites across " << _baselineElitesBySpecies.size() << " species" << std::endl;
    
    // Log initial elite state for debugging
    std::cout << "Initial elite distribution:" << std::endl;
    for (const auto& [speciesId, eliteIndices] : _baselineElitesBySpecies) {
        std::cout << "  Species " << speciesId << ": " << eliteIndices.size() << " elites" << std::endl;
        for (size_t eliteIdx : eliteIndices) {
            auto it = std::find(_baselineElites.begin(), _baselineElites.end(), eliteIdx);
            if (it != _baselineElites.end()) {
                size_t baselineIndex = std::distance(_baselineElites.begin(), it);
                std::cout << "    Elite " << eliteIdx << " fitness: " << _baselineEliteFitness[baselineIndex].getFitness() << std::endl;
            }
        }
    }
    
    // ========================================================================
    // PHASE 2: MULTI-GENERATION EVOLUTION - Test Elite Preservation
    // ========================================================================
    
    std::cout << "\n=== PHASE 2: MULTI-GENERATION EVOLUTION ===" << std::endl;
    
    const uint32_t NUM_TEST_GENERATIONS = 5;
    const uint32_t startingGeneration = _foundationGeneration + 1;
    
    for (uint32_t gen = startingGeneration; gen < startingGeneration + NUM_TEST_GENERATIONS; ++gen) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "TESTING GENERATION " << gen << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
        // Run a single generation of evolution with elite protection
        runSingleGenerationEvolution(gen);
        
        // ========================================================================
        // CRITICAL VALIDATIONS - Test the Exact Mechanisms Causing the Bug
        // ========================================================================
        
        // Validation 1: Elite Identity Preservation (EvolutionPrototype.hpp:441-448)
        validateEliteIdentityPreservation(gen);
        
        // Validation 2: Elite Fitness Preservation (EvolutionPrototype.hpp:544-565) 
        validateEliteFitnessPreservation(gen);
        
        // Validation 3: Elite Status Retention (GlobalIndexRegistry consistency)
        validateEliteStatusRetention(gen);
        
        // Validation 4: Species Fitness Progression (THE FAILING ASSERTION!)
        // This tests EvolutionPrototype.hpp:731-732 - the exact assertion that fails
        validateSpeciesFitnessProgression(gen);
        
        std::cout << "Generation " << gen << " validation completed successfully" << std::endl;
    }
    
    // ========================================================================
    // PHASE 3: COMPREHENSIVE VALIDATION - Final System State
    // ========================================================================
    
    std::cout << "\n=== PHASE 3: COMPREHENSIVE VALIDATION ===" << std::endl;
    
    const uint32_t finalGeneration = startingGeneration + NUM_TEST_GENERATIONS - 1;
    auto& finalGenomes = _populationContainer->getCurrentGenomes(finalGeneration);
    auto& finalFitnessResults = _populationContainer->getCurrentFitnessResults(finalGeneration);
    
    // Validate that all baseline elites still exist and are preserved
    std::cout << "Final elite preservation summary:" << std::endl;
    
    size_t preservedElites = 0;
    size_t identicalElites = 0;
    size_t fitnessPreservedElites = 0;
    
    for (size_t i = 0; i < _baselineElites.size(); ++i) {
        size_t eliteIndex = _baselineElites[i];
        
        // Check elite status preservation
        auto currentState = _globalIndexRegistry->getState(eliteIndex);
        bool stillElite = (currentState == GenomeState::Elite);
        if (stillElite) preservedElites++;
        
        // Check identity preservation
        GenomeProperties currentProps = extractGenomeProperties(finalGenomes[eliteIndex]);
        bool identityPreserved = currentProps.isIdenticalTo(_baselineEliteProperties[i]);
        if (identityPreserved) identicalElites++;
        
        // Check fitness preservation
        TestFitnessResult currentFitness(0.0);
        for (const auto& [fitness, globalIndex] : finalFitnessResults) {
            if (globalIndex == eliteIndex) {
                currentFitness = fitness;
                break;
            }
        }
        bool fitnessPreserved = currentFitness.isEqualTo(_baselineEliteFitness[i]);
        if (fitnessPreserved) fitnessPreservedElites++;
        
        std::cout << "  Elite " << eliteIndex << ": status=" << (stillElite ? "PRESERVED" : "LOST")
                  << ", identity=" << (identityPreserved ? "IDENTICAL" : "CHANGED")
                  << ", fitness=" << (fitnessPreserved ? "PRESERVED" : "CHANGED")
                  << " (" << _baselineEliteFitness[i].getFitness() << " -> " << currentFitness.getFitness() << ")"
                  << std::endl;
    }
    
    // Final success criteria
    EXPECT_EQ(preservedElites, _baselineElites.size()) 
        << "All elite genomes should retain Elite status across " << NUM_TEST_GENERATIONS << " generations";
    EXPECT_EQ(identicalElites, _baselineElites.size()) 
        << "All elite genomes should remain byte-for-byte identical across " << NUM_TEST_GENERATIONS << " generations";
    EXPECT_EQ(fitnessPreservedElites, _baselineElites.size()) 
        << "All elite fitness values should be preserved across " << NUM_TEST_GENERATIONS << " generations";
    
    // Species fitness progression validation (the critical assertion should never fail)
    std::cout << "\nSpecies fitness progression summary:" << std::endl;
    for (const auto& [speciesId, eliteIndices] : _baselineElitesBySpecies) {
        std::cout << "  Species " << speciesId << ": FITNESS PROGRESSION MAINTAINED (assertion never failed)" << std::endl;
    }
    
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "SUCCESS: Elite preservation mechanisms work correctly!" << std::endl;
    std::cout << "- Elite genomes preserved identically across " << NUM_TEST_GENERATIONS << " generations" << std::endl;
    std::cout << "- Elite fitness values maintained exactly" << std::endl;
    std::cout << "- Species fitness progression assertion never failed" << std::endl;
    std::cout << "- Elite status retained in GlobalIndexRegistry" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    
    // If this test passes, the species retention bug is NOT in the elite preservation mechanisms
    // The bug must be elsewhere (species elimination, population control, etc.)
}