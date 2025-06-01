#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/WeightMutation.hpp"
#include "version3/data/Genome.hpp"

using namespace Operator;

class WeightMutationTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::srand(42); // Fixed seed for reproducibility
    }

    // Helper to create a basic genome with connections for testing
    Genome createBasicGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3, 4};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, 
            {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1, 2, 3};
        params._sourceNodeHistoryIDs = {1, 2, 3};
        params._targetNodeHistoryIDs = {2, 3, 4};
        params._connectionAttributes = {
            {1.0f, true},   // enabled
            {-0.5f, true},  // enabled
            {2.0f, false}   // disabled
        };
        return Genome(params);
    }

    // Helper to create genome with only enabled connections
    Genome createEnabledOnlyGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1, 2};
        params._sourceNodeHistoryIDs = {1, 2};
        params._targetNodeHistoryIDs = {2, 3};
        params._connectionAttributes = {
            {1.5f, true},
            {-2.0f, true}
        };
        return Genome(params);
    }

    // Helper to create genome with no connections
    Genome createEmptyGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
        // No connections
        return Genome(params);
    }

    // Helper to create genome with all disabled connections
    Genome createAllDisabledGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1, 2};
        params._sourceNodeHistoryIDs = {1, 2};
        params._targetNodeHistoryIDs = {2, 3};
        params._connectionAttributes = {
            {1.0f, false},  // disabled
            {-1.0f, false}  // disabled
        };
        return Genome(params);
    }

    // Helper to count number of different weights between two genomes
    size_t countWeightDifferences(const Genome& genome1, const Genome& genome2) {
        const auto& conns1 = genome1.get_connectionGenes();
        const auto& conns2 = genome2.get_connectionGenes();
        
        if (conns1.size() != conns2.size()) return SIZE_MAX; // Invalid comparison
        
        size_t differences = 0;
        for (size_t i = 0; i < conns1.size(); ++i) {
            if (conns1[i].get_attributes().weight != conns2[i].get_attributes().weight) {
                differences++;
            }
        }
        return differences;
    }

    // Helper to verify genome structure integrity
    void verifyGenomeIntegrity(const Genome& original, const Genome& mutated) {
        // Same number of nodes and connections
        EXPECT_EQ(original.get_nodeGenes().size(), mutated.get_nodeGenes().size());
        EXPECT_EQ(original.get_connectionGenes().size(), mutated.get_connectionGenes().size());
        
        // Same node data
        for (size_t i = 0; i < original.get_nodeGenes().size(); ++i) {
            const auto& origNode = original.get_nodeGenes()[i];
            const auto& mutNode = mutated.get_nodeGenes()[i];
            EXPECT_EQ(origNode.get_historyID(), mutNode.get_historyID());
            EXPECT_EQ(origNode.get_type(), mutNode.get_type());
            EXPECT_EQ(origNode.get_attributes().activationType, mutNode.get_attributes().activationType);
        }
        
        // Same connection structure (but potentially different weights)
        for (size_t i = 0; i < original.get_connectionGenes().size(); ++i) {
            const auto& origConn = original.get_connectionGenes()[i];
            const auto& mutConn = mutated.get_connectionGenes()[i];
            EXPECT_EQ(origConn.get_historyID(), mutConn.get_historyID());
            EXPECT_EQ(origConn.get_sourceNodeGene().get_historyID(), 
                     mutConn.get_sourceNodeGene().get_historyID());
            EXPECT_EQ(origConn.get_targetNodeGene().get_historyID(), 
                     mutConn.get_targetNodeGene().get_historyID());
            EXPECT_EQ(origConn.get_attributes().enabled, mutConn.get_attributes().enabled);
        }
    }
};

// ============================================================================
// Weight Modification Verification Tests
// ============================================================================

TEST_F(WeightMutationTest, WeightActuallyChanges_PerturbationOnly) {
    Genome original = createEnabledOnlyGenome();
    WeightMutationParams params(1.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    
    Genome mutated = weightMutation(original, params);
    
    // With rate 1.0, all enabled connections should be mutated
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_GT(differences, 0u);
}

TEST_F(WeightMutationTest, WeightActuallyChanges_ReplacementOnly) {
    Genome original = createEnabledOnlyGenome();
    WeightMutationParams params(0.0, 1.0, 0.5, 2.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    
    Genome mutated = weightMutation(original, params);
    
    // With rate 1.0, all enabled connections should be mutated
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_GT(differences, 0u);
}

TEST_F(WeightMutationTest, WeightActuallyChanges_Mixed) {
    Genome original = createEnabledOnlyGenome();
    WeightMutationParams params(0.5, 0.5, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    Genome mutated = weightMutation(original, params);
    
    // With rates 0.5 each, some connections should be mutated
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_GT(differences, 0u);
}

TEST_F(WeightMutationTest, NoChangeWhenRateZero_Perturbation) {
    Genome original = createEnabledOnlyGenome();
    WeightMutationParams params(0.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    
    Genome mutated = weightMutation(original, params);
    
    // With rate 0.0, no weights should change
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_EQ(differences, 0u);
}

TEST_F(WeightMutationTest, NoChangeWhenRateZero_Replacement) {
    Genome original = createEnabledOnlyGenome();
    WeightMutationParams params(0.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    
    Genome mutated = weightMutation(original, params);
    
    // With rate 0.0, no weights should change
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_EQ(differences, 0u);
}

TEST_F(WeightMutationTest, DisabledConnectionsUntouched) {
    Genome original = createBasicGenome();
    WeightMutationParams params(1.0, 1.0, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    Genome mutated = weightMutation(original, params);
    
    const auto& origConns = original.get_connectionGenes();
    const auto& mutConns = mutated.get_connectionGenes();
    
    // Find disabled connection (should have same weight)
    for (size_t i = 0; i < origConns.size(); ++i) {
        if (!origConns[i].get_attributes().enabled) {
            EXPECT_FLOAT_EQ(origConns[i].get_attributes().weight, 
                           mutConns[i].get_attributes().weight);
        }
    }
}

TEST_F(WeightMutationTest, CorrectMutationType_PerturbationVsReplacement) {
    Genome original = createEnabledOnlyGenome();
    
    // Test perturbation - changes should be smaller
    WeightMutationParams perturbParams(1.0, 0.0, 0.1, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    Genome perturbMutated = weightMutation(original, perturbParams);
    
    // Test replacement - changes should be larger and within range
    WeightMutationParams replaceParams(0.0, 1.0, 0.1, 1.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    Genome replaceMutated = weightMutation(original, replaceParams);
    
    const auto& origConns = original.get_connectionGenes();
    const auto& perturbConns = perturbMutated.get_connectionGenes();
    const auto& replaceConns = replaceMutated.get_connectionGenes();
    
    for (size_t i = 0; i < origConns.size(); ++i) {
        if (origConns[i].get_attributes().enabled) {
            // Perturbation should be small relative to original
            float perturbDiff = std::abs(perturbConns[i].get_attributes().weight - 
                                       origConns[i].get_attributes().weight);
            EXPECT_LT(perturbDiff, 1.0f); // Should be small with strength 0.1
            
            // Replacement should be within specified range
            float replaceWeight = replaceConns[i].get_attributes().weight;
            EXPECT_GE(replaceWeight, -1.0f);
            EXPECT_LE(replaceWeight, 1.0f);
        }
    }
}

// ============================================================================
// Genome Integrity After Mutation Tests
// ============================================================================

TEST_F(WeightMutationTest, OriginalGenomeUntouched) {
    Genome original = createBasicGenome();
    
    // Store original weights
    std::vector<float> originalWeights;
    for (const auto& conn : original.get_connectionGenes()) {
        originalWeights.push_back(conn.get_attributes().weight);
    }
    
    WeightMutationParams params(1.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    Genome mutated = weightMutation(original, params);
    
    // Verify original weights unchanged
    const auto& origConns = original.get_connectionGenes();
    for (size_t i = 0; i < origConns.size(); ++i) {
        EXPECT_FLOAT_EQ(origConns[i].get_attributes().weight, originalWeights[i]);
    }
}

TEST_F(WeightMutationTest, StructurePreservation) {
    Genome original = createBasicGenome();
    WeightMutationParams params(0.5, 0.5, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    Genome mutated = weightMutation(original, params);
    
    verifyGenomeIntegrity(original, mutated);
}

TEST_F(WeightMutationTest, ValidPhenotypeConstruction) {
    Genome original = createBasicGenome();
    WeightMutationParams params(1.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    
    Genome mutated = weightMutation(original, params);
    
    // Should be able to construct phenotype without throwing
    EXPECT_NO_THROW({
        mutated.constructPhenotype();
        auto phenotype = mutated.get_phenotype();
        EXPECT_NE(phenotype, nullptr);
    });
}

TEST_F(WeightMutationTest, NodeGenesUnchanged) {
    Genome original = createBasicGenome();
    WeightMutationParams params(1.0, 1.0, 1.0, 3.0, WeightMutationParams::MutationType::MIXED);
    
    Genome mutated = weightMutation(original, params);
    
    const auto& origNodes = original.get_nodeGenes();
    const auto& mutNodes = mutated.get_nodeGenes();
    
    EXPECT_EQ(origNodes.size(), mutNodes.size());
    for (size_t i = 0; i < origNodes.size(); ++i) {
        EXPECT_EQ(origNodes[i].get_historyID(), mutNodes[i].get_historyID());
        EXPECT_EQ(origNodes[i].get_type(), mutNodes[i].get_type());
        EXPECT_EQ(origNodes[i].get_attributes().activationType, 
                 mutNodes[i].get_attributes().activationType);
    }
}

// ============================================================================
// Parameter Validation Tests  
// ============================================================================

TEST_F(WeightMutationTest, InvalidRateRejection_PerturbationNegative) {
    EXPECT_DEATH({
        WeightMutationParams params(-0.1, 0.5, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidRateRejection_PerturbationTooHigh) {
    EXPECT_DEATH({
        WeightMutationParams params(1.1, 0.5, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidRateRejection_ReplacementNegative) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, -0.1, 0.5, 2.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidRateRejection_ReplacementTooHigh) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, 1.1, 0.5, 2.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidStrengthRejection_Zero) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, 0.5, 0.0, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidStrengthRejection_Negative) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, 0.5, -0.1, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidRangeRejection_Zero) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, 0.5, 0.5, 0.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, InvalidRangeRejection_Negative) {
    EXPECT_DEATH({
        WeightMutationParams params(0.5, 0.5, 0.5, -1.0, WeightMutationParams::MutationType::REPLACEMENT_ONLY);
    }, "");
}

TEST_F(WeightMutationTest, ValidBoundaryValues) {
    // Test boundary values that should work
    EXPECT_NO_THROW({
        WeightMutationParams params1(0.0, 0.0, 0.001, 0.001, WeightMutationParams::MutationType::MIXED);
        WeightMutationParams params2(1.0, 1.0, 10.0, 100.0, WeightMutationParams::MutationType::MIXED);
    });
}

// ============================================================================
// Edge Case Handling Tests
// ============================================================================

TEST_F(WeightMutationTest, EmptyConnections) {
    Genome original = createEmptyGenome();
    WeightMutationParams params(1.0, 1.0, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    // Should return safely without errors
    EXPECT_NO_THROW({
        Genome mutated = weightMutation(original, params);
        EXPECT_EQ(mutated.get_connectionGenes().size(), 0u);
    });
}

TEST_F(WeightMutationTest, AllDisabledConnections) {
    Genome original = createAllDisabledGenome();
    WeightMutationParams params(1.0, 1.0, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    Genome mutated = weightMutation(original, params);
    
    // Output should be identical to input when all connections disabled
    size_t differences = countWeightDifferences(original, mutated);
    EXPECT_EQ(differences, 0u);
}

TEST_F(WeightMutationTest, SingleConnection) {
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2};
    params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    params._connectionHistoryIDs = {1};
    params._sourceNodeHistoryIDs = {1};
    params._targetNodeHistoryIDs = {2};
    params._connectionAttributes = {{1.0f, true}};
    
    Genome original(params);
    WeightMutationParams mutParams(1.0, 0.0, 0.5, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    
    // Should work correctly with minimal genome
    EXPECT_NO_THROW({
        Genome mutated = weightMutation(original, mutParams);
        EXPECT_EQ(mutated.get_connectionGenes().size(), 1u);
        
        // Weight should have changed
        float origWeight = original.get_connectionGenes()[0].get_attributes().weight;
        float mutWeight = mutated.get_connectionGenes()[0].get_attributes().weight;
        EXPECT_NE(origWeight, mutWeight);
    });
}

TEST_F(WeightMutationTest, LargeGenome_Performance) {
    // Create genome with many connections
    GenomeParams params;
    const size_t numNodes = 100;
    const size_t numConnections = 500;
    
    for (size_t i = 1; i <= numNodes; ++i) {
        params._nodeHistoryIDs.push_back(static_cast<uint32_t>(i));
        if (i <= 10) {
            params._nodeTypes.push_back(NodeType::INPUT);
        } else if (i > 90) {
            params._nodeTypes.push_back(NodeType::OUTPUT);
        } else {
            params._nodeTypes.push_back(NodeType::HIDDEN);
        }
        params._nodeAttributes.push_back({ActivationType::SIGMOID});
    }
    
    // Add many connections
    for (size_t i = 0; i < numConnections; ++i) {
        params._connectionHistoryIDs.push_back(static_cast<uint32_t>(i + 1));
        params._sourceNodeHistoryIDs.push_back(static_cast<uint32_t>((i % 50) + 1)); // Input/Hidden nodes
        params._targetNodeHistoryIDs.push_back(static_cast<uint32_t>((i % 50) + 51)); // Hidden/Output nodes
        params._connectionAttributes.push_back({static_cast<float>(i * 0.1), true});
    }
    
    Genome original(params);
    WeightMutationParams mutParams(0.1, 0.1, 0.5, 2.0, WeightMutationParams::MutationType::MIXED);
    
    // Should handle large genome efficiently
    EXPECT_NO_THROW({
        Genome mutated = weightMutation(original, mutParams);
        verifyGenomeIntegrity(original, mutated);
    });
}

TEST_F(WeightMutationTest, ExtremeWeightValues) {
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2};
    params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    params._connectionHistoryIDs = {1};
    params._sourceNodeHistoryIDs = {1};
    params._targetNodeHistoryIDs = {2};
    params._connectionAttributes = {{1000000.0f, true}}; // Very large weight
    
    Genome original(params);
    WeightMutationParams mutParams(1.0, 0.0, 0.1, 2.0, WeightMutationParams::MutationType::PERTURBATION_ONLY);
    
    // Should handle extreme values gracefully
    EXPECT_NO_THROW({
        Genome mutated = weightMutation(original, mutParams);
        EXPECT_EQ(mutated.get_connectionGenes().size(), 1u);
    });
}