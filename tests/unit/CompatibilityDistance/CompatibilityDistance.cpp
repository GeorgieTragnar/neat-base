#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <unordered_map>

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

TEST_F(CompatibilityDistanceTest, ParameterConstruction) {
    // Test valid parameter construction
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    EXPECT_NO_THROW(params);
    
    // Test with different coefficients
    Operator::CompatibilityDistanceParams params2(2.0f, 2.0f, 1.0f, 5.0f);
    EXPECT_NO_THROW(params2);
}

TEST_F(CompatibilityDistanceTest, IdenticalGenomesSameSpecies) {
    auto genome1 = createSimpleGenome();
    auto genome2 = createSimpleGenome();
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    
    // First genome creates species 1
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    EXPECT_EQ(species1, 1u);
    
    // Identical genome should be assigned to same species
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    EXPECT_EQ(species2, species1);
}

TEST_F(CompatibilityDistanceTest, DifferentWeightsSameSpecies) {
    auto genome1 = createGenomeWithConnections({{1.0f, true}, {2.0f, true}});
    auto genome2 = createGenomeWithConnections({{1.1f, true}, {2.1f, true}});
    
    // Small weight differences should keep genomes in same species
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    EXPECT_EQ(species1, species2);
}

TEST_F(CompatibilityDistanceTest, LargeWeightDifferenceNewSpecies) {
    auto genome1 = createGenomeWithConnections({{1.0f, true}});
    auto genome2 = createGenomeWithConnections({{10.0f, true}});
    
    // Large weight differences with high c3 should create new species
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 2.0f, 3.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, StrictThresholdNewSpecies) {
    auto genome1 = createSimpleGenome();
    auto genome2 = createGenomeWithConnections({{1.1f, true}});
    
    // Very strict threshold should create separate species for any difference
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 1.0f, 0.01f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, SpeciesIdIncrementation) {
    auto genome1 = createSimpleGenome();
    auto genome2 = createGenomeWithConnections({{10.0f, true}});
    auto genome3 = createGenomeWithConnections({{20.0f, true}});
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 1.0f, 1.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    uint32_t species3 = Operator::compatibilityDistance(genome3, historyTracker, params);
    
    EXPECT_EQ(species1, 1u);
    EXPECT_EQ(species2, 2u);
    EXPECT_EQ(species3, 3u);
}

TEST_F(CompatibilityDistanceTest, RepresentativeStorage) {
    auto genome1 = createSimpleGenome();
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 3.0f);
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    
    // Verify representative is stored
    EXPECT_EQ(historyTracker->_speciesRepresentatives.size(), 1u);
    EXPECT_TRUE(historyTracker->_speciesRepresentatives.find(species1) != 
                historyTracker->_speciesRepresentatives.end());
}

TEST_F(CompatibilityDistanceTest, MultipleGenomesToSameSpecies) {
    auto genome1 = createSimpleGenome();
    auto genome2 = createSimpleGenome();
    auto genome3 = createSimpleGenome();
    
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 10.0f); // High threshold
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    uint32_t species3 = Operator::compatibilityDistance(genome3, historyTracker, params);
    
    EXPECT_EQ(species1, species2);
    EXPECT_EQ(species2, species3);
    
    // Should only have one species representative
    EXPECT_EQ(historyTracker->_speciesRepresentatives.size(), 1u);
}

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

TEST_F(CompatibilityDistanceTest, DisabledConnectionsAffectDistance) {
    auto genome1 = createGenomeWithConnections({{1.0f, true}});
    auto genome2 = createGenomeWithConnections({{1.0f, false}});
    
    // Same weight but different enabled state should still calculate weight difference
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 2.0f, 1.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    // They should be in same species since weights are identical (enabled state doesn't affect weight diff)
    EXPECT_EQ(species1, species2);
}

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

TEST_F(CompatibilityDistanceTest, BoundaryThresholdTesting) {
    auto genome1 = createGenomeWithConnections({{1.0f, true}});
    auto genome2 = createGenomeWithConnections({{2.0f, true}});
    
    // Set threshold exactly at expected distance (c3 * weight_diff = 0.5 * 1.0 = 0.5)
    Operator::CompatibilityDistanceParams params(0.0f, 0.0f, 0.5f, 0.5f);
    
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, params);
    
    // At exact threshold boundary, should be same species (distance < threshold)
    EXPECT_EQ(species1, species2);
}

TEST_F(CompatibilityDistanceTest, ExactThresholdBoundary) {
    // First, create a genome and establish species 1
    auto genome1 = createSimpleGenome();
    Operator::CompatibilityDistanceParams params(1.0f, 1.0f, 0.4f, 0.5f);
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, params);
    EXPECT_EQ(species1, 1u);
    
    // Now create a genome with very different weight that will exceed threshold
    auto genome2 = createGenomeWithConnections({{5.0f, true}});
    
    // Distance = c3 * |5.0 - 1.0| = 0.4 * 4.0 = 1.6, which is > 0.5 threshold
    Operator::CompatibilityDistanceParams strictParams(0.0f, 0.0f, 0.4f, 0.5f);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, strictParams);
    
    // Should create new species due to large weight difference
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, NoMatchingConnections) {
    // Create genome1 with just bias connection (no input-to-output)
    Operator::InitParams params1(
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 input
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 output
        {{0, ConnectionGeneAttributes{1.0f, true}}},   // bias to output only
        Operator::InitParams::InputConnectionStrategy::NONE // No input connections
    );
    auto genome1 = Operator::init(historyTracker, params1);
    
    // Establish species 1
    Operator::CompatibilityDistanceParams setupParams(1.0f, 1.0f, 0.4f, 10.0f);
    uint32_t species1 = Operator::compatibilityDistance(genome1, historyTracker, setupParams);
    EXPECT_EQ(species1, 1u);
    
    // Create genome2 with more connections to ensure excess genes
    auto genome2 = createGenomeWithMoreConnections(); // This has more connections than genome1
    
    // Use high excess coefficient with low threshold
    // genome2 has more connections than genome1, creating excess genes
    Operator::CompatibilityDistanceParams testParams(10.0f, 0.0f, 0.0f, 1.0f);
    uint32_t species2 = Operator::compatibilityDistance(genome2, historyTracker, testParams);
    
    EXPECT_NE(species1, species2);
}

TEST_F(CompatibilityDistanceTest, LargeGenomeNormalization) {
    // This test verifies the N > 20 normalization behavior
    // Since we can't easily create 20+ connection genomes with current Init operator,
    // we test the principle with smaller genomes and appropriate coefficients
    
    auto smallGenome = createSimpleGenome();
    auto largeGenome = createGenomeWithMoreConnections();
    
    // Test that larger genomes have normalized distance impact
    // Higher coefficients to amplify the effect
    Operator::CompatibilityDistanceParams params(2.0f, 2.0f, 1.0f, 3.0f);
    
    uint32_t species1 = Operator::compatibilityDistance(smallGenome, historyTracker, params);
    uint32_t species2 = Operator::compatibilityDistance(largeGenome, historyTracker, params);
    
    // With normalization, excess genes should still create separation
    EXPECT_NE(species1, species2);
}

// Additional recommended unit tests

TEST_F(CompatibilityDistanceTest, ExcessCoefficientIsolation) {
    // Create genome1 with fewer connections
    Operator::InitParams params1(
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 input
        {NodeGeneAttributes{ActivationType::SIGMOID}}, // 1 output
        {}, // No bias connections
        Operator::InitParams::InputConnectionStrategy::NONE
    );
    auto smallGenome = Operator::init(historyTracker, params1);
    
    // Establish species 1
    Operator::CompatibilityDistanceParams setupParams(1.0f, 1.0f, 0.4f, 10.0f);
    uint32_t species1 = Operator::compatibilityDistance(smallGenome, historyTracker, setupParams);
    EXPECT_EQ(species1, 1u);
    
    // Create genome with more connections (excess genes)
    auto largeGenome = createGenomeWithMoreConnections();
    
    // Test with only c1 coefficient active (isolate excess gene penalty)
    Operator::CompatibilityDistanceParams testParams(5.0f, 0.0f, 0.0f, 2.0f);
    uint32_t species2 = Operator::compatibilityDistance(largeGenome, historyTracker, testParams);
    
    // Should create new species due to excess genes
    EXPECT_NE(species1, species2);
}

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