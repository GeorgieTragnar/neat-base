#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <unordered_map>
#include <iostream>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/CompatibilityDistance.hpp"
#include "version3/operator/Init.hpp"
#include "version3/data/HistoryTracker.hpp"

class CompatibilityDistanceTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_shared<HistoryTracker>();
    }
    
    Genome createSimpleGenome() {
        Operator::InitParams params(
            {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 input
            {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 output
            {{0, ConnectionGeneAttributes{1.0f, true}}},   // bias to output
            Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS
        );
        return Operator::init(historyTracker, params);
    }
    
    Genome createGenomeWithConnections(const std::vector<std::pair<float, bool>>& connectionData) {
        auto genome = createSimpleGenome();
        
        // Modify connections to have specific weights and enabled states
        auto& connections = genome.get_connectionGenes();
        for (size_t i = 0; i < std::min(connections.size(), connectionData.size()); ++i) {
            connections[i].get_attributes().weight = connectionData[i].first;
            connections[i].get_attributes().enabled = connectionData[i].second;
        }
        
        return genome;
    }
    
    Genome createGenomeWithMoreConnections() {
        // Create genome with multiple connections to test excess genes
        Operator::InitParams params(
            {NodeGeneAttributes{ActivationType::SIGMOID}, NodeGeneAttributes{ActivationType::SIGMOID}}, // 2 inputs
            {NodeGeneAttributes{ActivationType::SIGMOID}, NodeGeneAttributes{ActivationType::SIGMOID}}, // 2 outputs
            {{0, ConnectionGeneAttributes{1.0f, true}}, {1, ConnectionGeneAttributes{1.5f, true}}}, // bias connections
            Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS
        );
        return Operator::init(historyTracker, params);
    }

    std::shared_ptr<HistoryTracker> historyTracker;
};

// REMOVED: Redundant with ParameterValidation test below

// REMOVED: Redundant with EliteCopyPropagation test

// REMOVED: Redundant with ThresholdBehaviorInvariant test

// REMOVED: Redundant with ThresholdBehaviorInvariant test

// REMOVED: Redundant with ThresholdBehaviorInvariant test

// REMOVED: Redundant with SpeciesProgression test

// REMOVED: Redundant with SpeciesProgression test

// REMOVED: This test is redundant with IdenticalGenomesSameSpecies

// Critical missing tests

TEST_F(CompatibilityDistanceTest, ExcessGenesCreateNewSpecies) {
    auto smallGenome = createSimpleGenome();
    auto largeGenome = createGenomeWithMoreConnections();
    
    // High c1 coefficient should separate genomes with different connection counts
    Operator::CompatibilityDistanceParams params(3.0f, 1.0f, 0.1f, 2.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(smallGenome, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(largeGenome, historyTracker, params);
    
    EXPECT_NE(species1, species2);
}

// REMOVED: This test verifies implementation detail, not core invariant

TEST_F(CompatibilityDistanceTest, EmptyGenomeEdgeCase) {
    // Create minimal genome structure 
    Operator::InitParams params(
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 input
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 output
        {}, // No bias connections
        Operator::InitParams::InputConnectionStrategy::NONE // No input connections
    );
    auto emptyGenome = Operator::init(historyTracker, params);
    auto normalGenome = createSimpleGenome();
    
    // Empty vs normal genome should create separate species with high c1
    Operator::CompatibilityDistanceParams testParams(3.0f, 1.0f, 0.1f, 1.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(emptyGenome, historyTracker, testParams);
    uint32_t species2 = Operator::compatibilityDistance(normalGenome, historyTracker, testParams);
    
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, ZeroCoefficientTesting) {
    auto genome1 = createGenomeWithConnections({{1.0f, true}});
    auto genome2 = createGenomeWithConnections({{10.0f, true}});
    
    // Zero weight coefficient should ignore weight differences
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.0f, 2.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    // Should be same species since weight differences are ignored
    EXPECT_EQ(species1, species2);
}

// REMOVED: This test is redundant with ExactThresholdBoundary

TEST_F(CompatibilityDistanceTest, ExactThresholdBoundary) {
    // Test precise threshold boundary behavior with consistent parameters
    auto genome1 = createSimpleGenome();
    auto genome2 = createGenomeWithConnections({{5.0f, true}});
    
    // Use fresh historyTracker with consistent parameters
    // Distance = c3 * |5.0 - 1.0| = 0.4 * 4.0 = 1.6, which is > 0.5 threshold
    auto freshHistoryTracker = std::make_shared<HistoryTracker>();
    Operator::CompatibilityDistanceParams strictParams(0.0f, 0.0f, 0.4f, 0.5f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, freshHistoryTracker, strictParams);
    uint32_t species2 = Operator::compatibilityDistance(genome2, freshHistoryTracker, strictParams);
    
    // Should create new species due to large weight difference
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, NoMatchingConnections) {
    // Create genome1 with just bias connection (no input-to-output)
    auto freshHistoryTracker = std::make_shared<HistoryTracker>();
    Operator::InitParams params1(
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 input
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 output
        {{0, ConnectionGeneAttributes{1.0f, true}}},   // bias to output only
        Operator::InitParams::InputConnectionStrategy::NONE // No input connections
    );
    auto genome1 = Operator::init(freshHistoryTracker, params1);
    
    // Create genome2 with more connections to ensure excess genes
    auto genome2 = createGenomeWithMoreConnections(); // This has more connections than genome1
    
    // Use consistent parameters with high excess coefficient and low threshold
    // genome2 has more connections than genome1, creating excess genes
    Operator::CompatibilityDistanceParams testParams(10.0f, 0.0f, 0.0f, 1.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, freshHistoryTracker, testParams);
    uint32_t species2 = Operator::compatibilityDistance(genome2, freshHistoryTracker, testParams);
    
    EXPECT_NE(species1, species2);
}

// REMOVED: Can't properly test N>20 normalization with current Init operator

// Additional recommended unit tests

// REMOVED: This test is redundant with ExcessGenesCreateNewSpecies

TEST_F(CompatibilityDistanceTest, ZeroThreshold) {
    auto genome1 = createSimpleGenome();
    auto genome2 = createGenomeWithConnections({{1.1f, true}});
    
    // Zero threshold should create new species for any difference
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 1.0f, 0.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    // Any non-zero distance should exceed zero threshold
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, ParameterValidation) {
    // Test that parameter construction works with edge values
    EXPECT_NO_THROW(Operator::CompatibilityDistanceParams(0.0f, 0.0f, 0.0f, 0.0f));
    EXPECT_NO_THROW(Operator::CompatibilityDistanceParams(100.0f, 100.0f, 100.0f, 100.0f));
    
    // Test very small positive values
    EXPECT_NO_THROW(Operator::CompatibilityDistanceParams(0.001f, 0.001f, 0.001f, 0.001f));
}

// Core Invariant Tests

TEST_F(CompatibilityDistanceTest, DeterministicAssignment) {
    // Core invariant: same genome always gets same species assignment
    auto genome = createSimpleGenome();
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    
    // First assignment
    uint32_t species1 = Operator::compatibilityDistance(genome, historyTracker, params);
    
    // Reset history tracker and repeat - should get same result
    historyTracker = std::make_shared<HistoryTracker>();
    uint32_t species2 = Operator::compatibilityDistance(genome, historyTracker, params);
    
    EXPECT_EQ(species1, species2) << "Same genome should always get same species assignment";
}

TEST_F(CompatibilityDistanceTest, EliteCopyPropagation) {
    // Test the elite copy scenario: identical genomes should maintain same species
    auto originalGenome = createSimpleGenome();
    auto eliteCopy = createSimpleGenome(); // Functionally identical to original
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    
    // Assign original genome to species
    uint32_t originalSpecies = Operator::compatibilityDistance(originalGenome, historyTracker, params);
    
    // Elite copy should be assigned to same species
    uint32_t eliteSpecies = Operator::compatibilityDistance(eliteCopy, historyTracker, params);
    
    EXPECT_EQ(originalSpecies, eliteSpecies) << "Elite copies should maintain same species assignment";
    
    // Multiple elite copies should all get same species
    auto eliteCopy2 = createSimpleGenome();
    uint32_t eliteSpecies2 = Operator::compatibilityDistance(eliteCopy2, historyTracker, params);
    
    EXPECT_EQ(originalSpecies, eliteSpecies2) << "Multiple elite copies should all get same species";
}

TEST_F(CompatibilityDistanceTest, DistanceFormulaAccuracy) {
    // Test that the NEAT distance formula works correctly
    // Fix: modify BOTH connections to ensure larger weight difference
    auto genome1 = createGenomeWithConnections({{1.0f, true}, {1.0f, true}});
    auto genome2 = createGenomeWithConnections({{4.0f, true}, {4.0f, true}});
    
    // Test case 1: Average weight diff = (|4-1| + |4-1|)/2 = 3.0, should exceed threshold 1.5
    auto historyTracker1 = std::make_shared<HistoryTracker>();
    Operator::CompatibilityDistanceParams strictParams(0.0f, 0.0f, 1.0f, 1.5f); // threshold 1.5
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker1, strictParams);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker1, strictParams);
    
    EXPECT_NE(species1, species2) << "Distance 3.0 should exceed threshold 1.5";
    
    // Test case 2: Distance 3.0 should be below threshold 4.0
    auto historyTracker2 = std::make_shared<HistoryTracker>();
    Operator::CompatibilityDistanceParams lenientParams(0.0f, 0.0f, 1.0f, 4.0f);
    
    uint32_t species3 = Operator::compatibilityDistance(genome1, historyTracker2, lenientParams);
    uint32_t species4 = Operator::compatibilityDistance(genome2, historyTracker2, lenientParams);
    
    EXPECT_EQ(species3, species4) << "Distance 3.0 should be below threshold 4.0";
}

TEST_F(CompatibilityDistanceTest, ThresholdBehaviorInvariant) {
    // Test core threshold behavior: below threshold joins species, above creates new
    auto genome1 = createGenomeWithConnections({{1.0f, true}});
    auto genome2 = createGenomeWithConnections({{1.1f, true}});
    auto genome3 = createGenomeWithConnections({{10.0f, true}});
    
    Operator::CompatibilityDistanceParams params(0.0f, 0.0f, 1.0f, 0.5f); // weight diff threshold
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    
    // Small weight difference: 1.0 * |1.1 - 1.0| = 0.1 < 0.5 threshold
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    EXPECT_EQ(species1, species2) << "Below threshold should join existing species";
    
    // Large weight difference: 1.0 * |10.0 - 1.0| = 9.0 > 0.5 threshold
    uint32_t species3 = Operator::compatibilityDistance(genome3, historyTracker, params);
    EXPECT_NE(species1, species3) << "Above threshold should create new species";
}

TEST_F(CompatibilityDistanceTest, SpeciesProgression) {
    // Test that species IDs progress sequentially and representatives are stored
    auto genome1 = createSimpleGenome();
    auto genome2 = createGenomeWithConnections({{10.0f, true}});
    auto genome3 = createGenomeWithConnections({{20.0f, true}});
    
    // Strict parameters to ensure each genome gets its own species
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 1.0f, 0.1f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    uint32_t species3 = Operator::compatibilityDistance(genome3, historyTracker, params);
    
    // Species IDs should progress sequentially
    EXPECT_EQ(species1, 0u);
    EXPECT_EQ(species2, 1u);
    EXPECT_EQ(species3, 2u);
    
    // Should have 3 representatives stored
    EXPECT_EQ(historyTracker->_speciesRepresentatives.size(), 3u);
}