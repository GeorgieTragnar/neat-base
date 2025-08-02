#include <gtest/gtest.h>
#include "version3/data/Genome.hpp"
#include "version3/data/PopulationData.hpp"
#include "version3/data/PopulationContainer.hpp"

using namespace Population;

// Simple fitness result type for testing
struct TestFitness {
    double value;
    bool operator<(const TestFitness& other) const { return value < other.value; }
    TestFitness(double v) : value(v) {}
};

class PopulationContainerTest : public ::testing::Test {
protected:
    void SetUp() override {
        registry = std::make_unique<GlobalIndexRegistry>(0);
        container = std::make_unique<PopulationContainer<TestFitness>>(*registry);
    }
    
    // Helper to create a simple genome for testing
    Genome createTestGenome(uint32_t inputNodeId = 1, uint32_t outputNodeId = 2, uint32_t connHistoryId = 10) {
        GenomeParams params;
        params._nodeHistoryIDs = {inputNodeId, outputNodeId};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
        params._connectionHistoryIDs = {connHistoryId};
        params._sourceNodeHistoryIDs = {inputNodeId};
        params._targetNodeHistoryIDs = {outputNodeId};
        params._connectionAttributes = {{1.0f, true}};
        return Genome(params);
    }
    
    // Helper to create test genome data
    DynamicGenomeData createTestGenomeData(uint32_t speciesId = 1) {
        DynamicGenomeData data;
        data.speciesId = speciesId;
        data.pendingEliminationCounter = 0;
        data.repairAttempts = 0;
        data.isUnderRepair = false;
        data.isMarkedForElimination = false;
        data.genomeIndex = UINT32_MAX;
        data.parentAIndex = UINT32_MAX;
        data.parentBIndex = UINT32_MAX;
        return data;
    }
    
    std::unique_ptr<GlobalIndexRegistry> registry;
    std::unique_ptr<PopulationContainer<TestFitness>> container;
};

// ============================================================================
// Generation Access and Modulo Mapping Tests
// ============================================================================

TEST_F(PopulationContainerTest, GenerationAccess_ModuloMapping) {
    // Test: Generation access correctly maps to triple-buffer indices
    
    // Test basic modulo mapping for genomes
    auto& gen0 = container->getGenomes(0);
    auto& gen1 = container->getGenomes(1);
    auto& gen2 = container->getGenomes(2);
    auto& gen3 = container->getGenomes(3);  // Should map to buffer 0
    auto& gen4 = container->getGenomes(4);  // Should map to buffer 1
    auto& gen5 = container->getGenomes(5);  // Should map to buffer 2
    
    // Verify same buffer for generation % 3
    EXPECT_EQ(&gen0, &gen3);
    EXPECT_EQ(&gen1, &gen4);
    EXPECT_EQ(&gen2, &gen5);
    
    // Test for genome data
    auto& data0 = container->getGenomeData(0);
    auto& data3 = container->getGenomeData(3);
    EXPECT_EQ(&data0, &data3);
    
    // Test for fitness results
    auto& fitness0 = container->getFitnessResults(0);
    auto& fitness3 = container->getFitnessResults(3);
    EXPECT_EQ(&fitness0, &fitness3);
}

TEST_F(PopulationContainerTest, GenerationAccess_ConsistentMapping) {
    // Test: Const and non-const accessors return same underlying data
    
    const auto& constContainer = *container;
    
    // Compare const and non-const genome access
    auto& genomes = container->getGenomes(5);
    const auto& constGenomes = constContainer.getGenomes(5);
    EXPECT_EQ(&genomes, &constGenomes);
    
    // Compare const and non-const genome data access
    auto& genomeData = container->getGenomeData(5);
    const auto& constGenomeData = constContainer.getGenomeData(5);
    EXPECT_EQ(&genomeData, &constGenomeData);
    
    // Compare const and non-const fitness results access
    auto& fitnessResults = container->getFitnessResults(5);
    const auto& constFitnessResults = constContainer.getFitnessResults(5);
    EXPECT_EQ(&fitnessResults, &constFitnessResults);
}

// ============================================================================
// Underflow Protection Tests
// ============================================================================

TEST_F(PopulationContainerTest, UnderflowProtection_AccessorAssertions) {
    // Test: Last generation accessors properly assert on underflow
    
    // Generation 0 - cannot access last generation
    EXPECT_DEATH(container->getLastGenomes(0), "Cannot access last generation: currentGen underflow");
    EXPECT_DEATH(container->getLastGenomeData(0), "Cannot access last generation data: currentGen underflow");
    EXPECT_DEATH(container->getLastFitnessResults(0), "Cannot access last generation fitness: currentGen underflow");
    
    // Generation 1 - cannot access generation before last
    EXPECT_DEATH(container->getGenerationBeforeLastGenomes(1), "Cannot access generation before last: currentGen underflow");
    EXPECT_DEATH(container->getGenerationBeforeLastGenomeData(1), "Cannot access generation before last data: currentGen underflow");
    EXPECT_DEATH(container->getGenerationBeforeLastFitnessResults(1), "Cannot access generation before last fitness: currentGen underflow");
}

TEST_F(PopulationContainerTest, UnderflowProtection_ValidAccess) {
    // Test: Valid generation access works correctly
    
    // Generation 1 - can access last (gen 0)
    auto& lastGenomes = container->getLastGenomes(1);
    auto& gen0Genomes = container->getGenomes(0);
    EXPECT_EQ(&lastGenomes, &gen0Genomes);
    
    // Generation 2 - can access last (gen 1) and before last (gen 0)
    auto& gen2LastGenomes = container->getLastGenomes(2);
    auto& gen1Genomes = container->getGenomes(1);
    EXPECT_EQ(&gen2LastGenomes, &gen1Genomes);
    
    auto& beforeLastGenomes = container->getGenerationBeforeLastGenomes(2);
    EXPECT_EQ(&beforeLastGenomes, &gen0Genomes);
}

TEST_F(PopulationContainerTest, UnderflowProtection_HighGenerations) {
    // Test: High generation numbers with underflow protection
    
    uint32_t currentGen = 100;
    
    // Should work fine with high generation numbers
    auto& current = container->getCurrentGenomes(currentGen);
    auto& last = container->getLastGenomes(currentGen);
    auto& beforeLast = container->getGenerationBeforeLastGenomes(currentGen);
    
    // Verify correct modulo mapping
    auto& expected100 = container->getGenomes(100);  // 100 % 3 = 1
    auto& expected99 = container->getGenomes(99);    // 99 % 3 = 0  
    auto& expected98 = container->getGenomes(98);    // 98 % 3 = 2
    
    EXPECT_EQ(&current, &expected100);
    EXPECT_EQ(&last, &expected99);
    EXPECT_EQ(&beforeLast, &expected98);
}

// ============================================================================
// Size Synchronization Tests
// ============================================================================

TEST_F(PopulationContainerTest, SizeSynchronization_InitialState) {
    // Test: Initial state has consistent empty sizes
    
    EXPECT_TRUE(container->validateConsistency());
    EXPECT_EQ(container->getGenerationSize(0), 0u);
    EXPECT_EQ(container->getGenerationSize(1), 0u);
    EXPECT_EQ(container->getGenerationSize(2), 0u);
}

TEST_F(PopulationContainerTest, SizeSynchronization_AfterPushBack) {
    // Test: push_back maintains size consistency across all generations
    
    auto genome = createTestGenome();
    auto genomeData = createTestGenomeData();
    
    // Add to generation 0
    uint32_t globalIndex = container->push_back(0, std::move(genome), std::move(genomeData));
    EXPECT_EQ(globalIndex, 0u);  // First genome gets index 0
    
    // All generations should have same size
    EXPECT_TRUE(container->validateConsistency());
    EXPECT_EQ(container->getGenerationSize(0), 1u);
    EXPECT_EQ(container->getGenerationSize(1), 1u);
    EXPECT_EQ(container->getGenerationSize(2), 1u);
    
    // Add to generation 1
    auto genome2 = createTestGenome(3, 4, 11);
    auto genomeData2 = createTestGenomeData(2);
    uint32_t globalIndex2 = container->push_back(1, std::move(genome2), std::move(genomeData2));
    EXPECT_EQ(globalIndex2, 1u);  // Second genome gets index 1
    
    // All generations should still have same size
    EXPECT_TRUE(container->validateConsistency());
    EXPECT_EQ(container->getGenerationSize(0), 2u);
    EXPECT_EQ(container->getGenerationSize(1), 2u);
    EXPECT_EQ(container->getGenerationSize(2), 2u);
}

TEST_F(PopulationContainerTest, SizeSynchronization_CrossGenerationConsistency) {
    // Test: Cross-generation consistency after multiple additions
    
    // Add genomes to different generations
    for (int i = 0; i < 5; ++i) {
        uint32_t targetGeneration = i % 3;  // Distribute across generations
        auto genome = createTestGenome(i * 2 + 1, i * 2 + 2, i + 10);
        auto genomeData = createTestGenomeData(i + 1);
        
        uint32_t globalIndex = container->push_back(targetGeneration, std::move(genome), std::move(genomeData));
        EXPECT_EQ(globalIndex, static_cast<uint32_t>(i));
        
        // Consistency should be maintained after each addition
        EXPECT_TRUE(container->validateConsistency());
        
        size_t expectedSize = i + 1;
        EXPECT_EQ(container->getGenerationSize(0), expectedSize);
        EXPECT_EQ(container->getGenerationSize(1), expectedSize);
        EXPECT_EQ(container->getGenerationSize(2), expectedSize);
    }
}

// ============================================================================
// Container Operations Integrity Tests
// ============================================================================

TEST_F(PopulationContainerTest, PushBackIntegrity_SynchronizedAddition) {
    // Test: push_back properly synchronizes with registry and maintains vector consistency
    
    auto genome = createTestGenome();
    auto genomeData = createTestGenomeData();
    
    // Registry should start with size 0
    EXPECT_EQ(registry->getMaxIndex(), 0u);
    
    // Push back should increment registry and return correct index
    uint32_t globalIndex = container->push_back(0, std::move(genome), std::move(genomeData));
    EXPECT_EQ(globalIndex, 0u);  // Should be the first index
    EXPECT_EQ(registry->getMaxIndex(), 1u);  // Registry should grow
    
    // Container should have consistent sizes
    EXPECT_TRUE(container->validateConsistency());
    EXPECT_EQ(container->getGenerationSize(0), 1u);
    EXPECT_EQ(container->getGenerationSize(1), 1u);
    EXPECT_EQ(container->getGenerationSize(2), 1u);
}

TEST_F(PopulationContainerTest, PushBackIntegrity_SequentialIndices) {
    // Test: Sequential push_back operations return sequential indices
    
    std::vector<uint32_t> indices;
    uint32_t initialRegistrySize = registry->getMaxIndex();
    
    // Add multiple genomes
    for (int i = 0; i < 3; ++i) {
        auto genome = createTestGenome(i * 2 + 1, i * 2 + 2, i + 10);
        auto genomeData = createTestGenomeData(i + 1);
        
        uint32_t globalIndex = container->push_back(i % 3, std::move(genome), std::move(genomeData));
        indices.push_back(globalIndex);
    }
    
    // Indices should be sequential starting from initial registry size
    for (size_t i = 0; i < indices.size(); ++i) {
        EXPECT_EQ(indices[i], initialRegistrySize + i);
    }
    
    // Registry should have grown appropriately
    EXPECT_EQ(registry->getMaxIndex(), initialRegistrySize + 3);
}

TEST_F(PopulationContainerTest, CapacityManagement_ReserveCapacity) {
    // Test: reserveCapacity properly reserves for all three buffers
    
    size_t reserveSize = 100;
    container->reserveCapacity(reserveSize);
    
    // All generations should have reserved capacity
    EXPECT_GE(container->getGenerationCapacity(0), reserveSize);
    EXPECT_GE(container->getGenerationCapacity(1), reserveSize);
    EXPECT_GE(container->getGenerationCapacity(2), reserveSize);
}

TEST_F(PopulationContainerTest, FitnessResultsManagement_ClearOperation) {
    // Test: clearGenerationFitnessResults works correctly
    
    // Add some fitness results
    auto& fitness0 = container->getFitnessResults(0);
    auto& fitness1 = container->getFitnessResults(1);
    fitness0.insert({TestFitness(10.0), 0});
    fitness1.insert({TestFitness(20.0), 1});
    
    EXPECT_EQ(fitness0.size(), 1u);
    EXPECT_EQ(fitness1.size(), 1u);
    
    // Clear generation 0
    container->clearGenerationFitnessResults(0);
    
    EXPECT_EQ(fitness0.size(), 0u);
    EXPECT_EQ(fitness1.size(), 1u);  // Generation 1 should be unaffected
}

// ============================================================================
// Validation Methods Tests
// ============================================================================

TEST_F(PopulationContainerTest, ValidationMethods_ConsistencyDetection) {
    // Test: validateConsistency accurately detects size mismatches
    
    // Initially consistent
    EXPECT_TRUE(container->validateConsistency());
    
    // Add genome normally - should remain consistent
    auto genome = createTestGenome();
    auto genomeData = createTestGenomeData();
    container->push_back(0, std::move(genome), std::move(genomeData));
    EXPECT_TRUE(container->validateConsistency());
    
    // Manually break consistency by directly accessing vectors
    // (This simulates a bug that validation should catch)
    auto& genomes = container->getGenomes(0);
    auto& genomeData0 = container->getGenomeData(0);
    
    // Force size mismatch
    genomes.pop_back();
    EXPECT_FALSE(container->validateConsistency());
    
    // Restore consistency
    auto genome2 = createTestGenome(3, 4, 11);
    genomes.emplace_back(std::move(genome2));
    EXPECT_TRUE(container->validateConsistency());
}

// ============================================================================
// Triple-Buffer Rotation Tests
// ============================================================================

TEST_F(PopulationContainerTest, TripleBufferRotation_TemporalAccess) {
    // Test: Temporal access patterns work correctly across generation progression
    
    // Simulate generation progression
    for (uint32_t currentGen = 0; currentGen < 10; ++currentGen) {
        // Add genome to current generation
        auto genome = createTestGenome(currentGen * 2 + 1, currentGen * 2 + 2, currentGen + 10);
        auto genomeData = createTestGenomeData(currentGen + 1);
        container->push_back(currentGen, std::move(genome), std::move(genomeData));
        
        // Test current generation access
        auto& current = container->getCurrentGenomes(currentGen);
        EXPECT_EQ(current.size(), currentGen + 1);
        
        // Test last generation access (if valid)
        if (currentGen > 0) {
            auto& last = container->getLastGenomes(currentGen);
            EXPECT_EQ(last.size(), currentGen + 1);  // Same size due to cross-generation sync
        }
        
        // Test generation before last access (if valid)
        if (currentGen > 1) {
            auto& beforeLast = container->getGenerationBeforeLastGenomes(currentGen);
            EXPECT_EQ(beforeLast.size(), currentGen + 1);  // Same size due to cross-generation sync
        }
    }
}

TEST_F(PopulationContainerTest, TripleBufferRotation_ModuloWrapping) {
    // Test: Buffer rotation works correctly when generation numbers wrap around modulo 3
    
    // Test critical transition points
    std::vector<uint32_t> testGenerations = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 1000};
    
    for (uint32_t gen : testGenerations) {
        auto& genomes = container->getGenomes(gen);
        auto& expectedBuffer = container->getGenomes(gen % 3);
        EXPECT_EQ(&genomes, &expectedBuffer) << "Failed at generation " << gen;
        
        // Test temporal access if valid
        if (gen > 1) {
            auto& last = container->getLastGenomes(gen);
            auto& expectedLast = container->getGenomes((gen - 1) % 3);
            EXPECT_EQ(&last, &expectedLast) << "Last generation failed at generation " << gen;
            
            auto& beforeLast = container->getGenerationBeforeLastGenomes(gen);
            auto& expectedBeforeLast = container->getGenomes((gen - 2) % 3);
            EXPECT_EQ(&beforeLast, &expectedBeforeLast) << "Before last generation failed at generation " << gen;
        }
    }
}