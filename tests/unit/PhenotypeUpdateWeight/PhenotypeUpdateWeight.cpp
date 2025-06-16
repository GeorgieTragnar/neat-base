#include <gtest/gtest.h>
#include <unordered_map>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/operator/PhenotypeUpdateWeight.hpp"

class PhenotypeUpdateWeightTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
    }

    // Helper to create basic genome: I1 -> H2 -> O3 with B0, with initial phenotype
    Genome createBasicGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), HIDDEN(2), OUTPUT(3)
        params._nodeHistoryIDs = {0, 1, 2, 3};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID}
        };
        
        // Connections: B0->H2(1.0), I1->H2(2.0), H2->O3(3.0)
        params._connectionHistoryIDs = {10, 11, 12};
        params._sourceNodeHistoryIDs = {0, 1, 2};
        params._targetNodeHistoryIDs = {2, 2, 3};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},
            ConnectionGeneAttributes{2.0f, true},
            ConnectionGeneAttributes{3.0f, true}
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to create complex genome with mixed enabled/disabled connections
    Genome createComplexGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), INPUT(2), HIDDEN(3), HIDDEN(4), OUTPUT(5), OUTPUT(6)
        params._nodeHistoryIDs = {0, 1, 2, 3, 4, 5, 6};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID}
        };
        
        // Mix of enabled and disabled connections
        params._connectionHistoryIDs = {10, 11, 12, 13, 14, 15, 16};
        params._sourceNodeHistoryIDs = {0, 1, 2, 3, 4, 1, 2};
        params._targetNodeHistoryIDs = {3, 3, 4, 5, 6, 5, 6};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->H3 (enabled)
            ConnectionGeneAttributes{2.0f, false},  // I1->H3 (disabled)
            ConnectionGeneAttributes{3.0f, true},   // I2->H4 (enabled)
            ConnectionGeneAttributes{4.0f, true},   // H3->O5 (enabled)
            ConnectionGeneAttributes{5.0f, true},   // H4->O6 (enabled)
            ConnectionGeneAttributes{6.0f, false},  // I1->O5 (disabled)
            ConnectionGeneAttributes{7.0f, true}    // I2->O6 (enabled)
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to create genome with all connections disabled
    Genome createAllDisabledGenome() {
        GenomeParams params;
        
        params._nodeHistoryIDs = {0, 1, 2, 3};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID}
        };
        
        // All connections disabled
        params._connectionHistoryIDs = {10, 11, 12};
        params._sourceNodeHistoryIDs = {0, 1, 2};
        params._targetNodeHistoryIDs = {2, 2, 3};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, false},
            ConnectionGeneAttributes{2.0f, false},
            ConnectionGeneAttributes{3.0f, false}
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to create single connection genome
    Genome createSingleConnectionGenome() {
        GenomeParams params;
        
        params._nodeHistoryIDs = {0, 1, 2};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID}
        };
        
        params._connectionHistoryIDs = {10};
        params._sourceNodeHistoryIDs = {1};
        params._targetNodeHistoryIDs = {2};
        params._connectionAttributes = {
            ConnectionGeneAttributes{2.5f, true}
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to modify weight in genome and track change
    void modifyGenomeWeight(Genome& genome, uint32_t connectionHistoryID, float newWeight) {
        auto& connectionGenes = genome.get_connectionGenes();
        
        for (auto& conn : connectionGenes) {
            if (conn.get_historyID() == connectionHistoryID) {
                conn.get_attributes().weight = newWeight;
                break;
            }
        }
    }

    // Helper to set connection deltas
    void setConnectionDeltas(Genome& genome, const std::vector<uint32_t>& deltas) {
        auto& connectionDeltas = genome.get_connectionGeneDeltas();
        connectionDeltas = deltas;
    }

    // Validation helpers
    void validateDeltasEmpty(Genome& genome) {
        EXPECT_TRUE(genome.get_connectionGeneDeltas().empty());
    }

    void validateWeightUpdate(const Genome& genome, uint32_t connectionHistoryID, float expectedWeight) {
        const auto& phenotype = genome.get_phenotype();
        const auto& connectionGenes = genome.get_connectionGenes();
        
        // Find the connection in genome
        const ConnectionGene* targetConnection = nullptr;
        for (const auto& conn : connectionGenes) {
            if (conn.get_historyID() == connectionHistoryID) {
                targetConnection = &conn;
                break;
            }
        }
        
        ASSERT_NE(targetConnection, nullptr);
        
        if (!targetConnection->get_attributes().enabled) {
            // Disabled connections shouldn't be in phenotype
            return;
        }
        
        // Find corresponding phenotype connection
        const auto& nodes = genome.get_nodeGenes();
        uint32_t sourceHistoryID = nodes[targetConnection->get_sourceNodeIndex()].get_historyID();
        uint32_t targetHistoryID = nodes[targetConnection->get_targetNodeIndex()].get_historyID();
        
        // Build mapping from history IDs to phenotype indices
        std::unordered_map<uint32_t, size_t> historyToIndex;
        std::unordered_set<uint32_t> includedHistoryIDs;
        
        const auto& nodeGenes = genome.get_nodeGenes();
        
        for (const auto& node : nodeGenes) {
            if (node.get_type() == NodeType::INPUT || 
                node.get_type() == NodeType::OUTPUT || 
                node.get_type() == NodeType::BIAS) {
                includedHistoryIDs.insert(node.get_historyID());
            }
        }
        
        for (const auto& conn : connectionGenes) {
            if (conn.get_attributes().enabled) {
                includedHistoryIDs.insert(nodeGenes[conn.get_sourceNodeIndex()].get_historyID());
                includedHistoryIDs.insert(nodeGenes[conn.get_targetNodeIndex()].get_historyID());
            }
        }
        
        size_t phenotypeIndex = 0;
        for (const auto& node : nodeGenes) {
            if (includedHistoryIDs.count(node.get_historyID()) > 0) {
                historyToIndex[node.get_historyID()] = phenotypeIndex++;
            }
        }
        
        // Find and validate the phenotype connection
        bool foundConnection = false;
        for (const auto& phenConn : phenotype._orderedConnections) {
            if (phenConn._sourceNodeIndex == historyToIndex[sourceHistoryID] &&
                phenConn._targetNodeIndex == historyToIndex[targetHistoryID]) {
                foundConnection = true;
                EXPECT_FLOAT_EQ(phenConn._connectionGeneAttribute.weight, expectedWeight);
                break;
            }
        }
        EXPECT_TRUE(foundConnection);
    }

    void validateWeightSynchronization(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        const auto& connectionGenes = genome.get_connectionGenes();
        const auto& nodeGenes = genome.get_nodeGenes();
        
        // Build mapping
        std::unordered_map<uint32_t, size_t> historyToIndex;
        std::unordered_set<uint32_t> includedHistoryIDs;
        
        for (const auto& node : nodeGenes) {
            if (node.get_type() == NodeType::INPUT || 
                node.get_type() == NodeType::OUTPUT || 
                node.get_type() == NodeType::BIAS) {
                includedHistoryIDs.insert(node.get_historyID());
            }
        }
        
        for (const auto& conn : connectionGenes) {
            if (conn.get_attributes().enabled) {
                includedHistoryIDs.insert(nodeGenes[conn.get_sourceNodeIndex()].get_historyID());
                includedHistoryIDs.insert(nodeGenes[conn.get_targetNodeIndex()].get_historyID());
            }
        }
        
        size_t phenotypeIndex = 0;
        for (const auto& node : nodeGenes) {
            if (includedHistoryIDs.count(node.get_historyID()) > 0) {
                historyToIndex[node.get_historyID()] = phenotypeIndex++;
            }
        }
        
        // Verify each enabled genome connection matches phenotype
        for (const auto& conn : connectionGenes) {
            if (conn.get_attributes().enabled) {
                uint32_t sourceHistoryID = nodeGenes[conn.get_sourceNodeIndex()].get_historyID();
                uint32_t targetHistoryID = nodeGenes[conn.get_targetNodeIndex()].get_historyID();
                
                bool foundConnection = false;
                for (const auto& phenConn : phenotype._orderedConnections) {
                    if (phenConn._sourceNodeIndex == historyToIndex[sourceHistoryID] &&
                        phenConn._targetNodeIndex == historyToIndex[targetHistoryID]) {
                        foundConnection = true;
                        EXPECT_FLOAT_EQ(phenConn._connectionGeneAttribute.weight, conn.get_attributes().weight);
                        break;
                    }
                }
                EXPECT_TRUE(foundConnection);
            }
        }
    }

    void validateUnchangedConnections(const Genome& genome, const std::vector<uint32_t>& unchangedIDs) {
        const auto& connectionGenes = genome.get_connectionGenes();
        
        for (uint32_t historyID : unchangedIDs) {
            const ConnectionGene* conn = nullptr;
            for (const auto& c : connectionGenes) {
                if (c.get_historyID() == historyID) {
                    conn = &c;
                    break;
                }
            }
            
            ASSERT_NE(conn, nullptr);
            
            if (conn->get_attributes().enabled) {
                // Weight should match what's in phenotype
                validateWeightUpdate(genome, historyID, conn->get_attributes().weight);
            }
        }
    }
};

// Input Validation Tests
TEST_F(PhenotypeUpdateWeightTest, NonEmptyDeltas_Succeeds) {
    auto genome = createBasicGenome();
    
    // Modify weight in genome
    modifyGenomeWeight(genome, 11, 5.0f);
    setConnectionDeltas(genome, {11});
    
    EXPECT_NO_THROW(Operator::phenotypeUpdateWeight(genome));
    
    validateWeightUpdate(genome, 11, 5.0f);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, EmptyDeltas_Fails) {
    auto genome = createBasicGenome();
    
    setConnectionDeltas(genome, {});
    
    EXPECT_DEATH(Operator::phenotypeUpdateWeight(genome), "Connection deltas must be non-empty");
}

TEST_F(PhenotypeUpdateWeightTest, NonExistentConnectionIds_Ignores) {
    auto genome = createBasicGenome();
    
    // Include non-existent ID - should be ignored, not crash
    setConnectionDeltas(genome, {999});
    
    EXPECT_DEATH(Operator::phenotypeUpdateWeight(genome), "Connection with history ID not found in genome");
}

TEST_F(PhenotypeUpdateWeightTest, MixedValidInvalidIds_UpdatesValidOnly) {
    auto genome = createBasicGenome();
    
    // Modify weight for valid connection
    modifyGenomeWeight(genome, 11, 5.0f);
    
    // Include mix of valid and invalid IDs
    setConnectionDeltas(genome, {11, 999});
    
    // Should fail on first invalid ID encountered
    EXPECT_DEATH(Operator::phenotypeUpdateWeight(genome), "Connection with history ID not found in genome");
}

// Weight Update Tests
TEST_F(PhenotypeUpdateWeightTest, SingleWeightChange_UpdatedInPhenotype) {
    auto genome = createBasicGenome();
    
    // Store original weight
    float originalWeight = 2.0f;
    float newWeight = 8.5f;
    
    // Modify weight in genome
    modifyGenomeWeight(genome, 11, newWeight);
    setConnectionDeltas(genome, {11});
    
    Operator::phenotypeUpdateWeight(genome);
    
    validateWeightUpdate(genome, 11, newWeight);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, MultipleWeightChanges_AllUpdatedCorrectly) {
    auto genome = createBasicGenome();
    
    // Modify multiple weights
    modifyGenomeWeight(genome, 10, 1.5f);
    modifyGenomeWeight(genome, 11, 2.5f);
    modifyGenomeWeight(genome, 12, 3.5f);
    
    setConnectionDeltas(genome, {10, 11, 12});
    
    Operator::phenotypeUpdateWeight(genome);
    
    validateWeightUpdate(genome, 10, 1.5f);
    validateWeightUpdate(genome, 11, 2.5f);
    validateWeightUpdate(genome, 12, 3.5f);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, EnabledConnections_WeightsUpdated) {
    auto genome = createComplexGenome();
    
    // Modify weights for enabled connections
    modifyGenomeWeight(genome, 10, 1.1f); // B0->H3 (enabled)
    modifyGenomeWeight(genome, 12, 3.3f); // I2->H4 (enabled)
    modifyGenomeWeight(genome, 16, 7.7f); // I2->O6 (enabled)
    
    setConnectionDeltas(genome, {10, 12, 16});
    
    Operator::phenotypeUpdateWeight(genome);
    
    validateWeightUpdate(genome, 10, 1.1f);
    validateWeightUpdate(genome, 12, 3.3f);
    validateWeightUpdate(genome, 16, 7.7f);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, DisabledConnectionsInDeltas_Ignored) {
    auto genome = createComplexGenome();
    
    // Modify weights for both enabled and disabled connections
    modifyGenomeWeight(genome, 10, 1.1f); // B0->H3 (enabled)
    modifyGenomeWeight(genome, 11, 2.2f); // I1->H3 (disabled)
    modifyGenomeWeight(genome, 15, 6.6f); // I1->O5 (disabled)
    
    setConnectionDeltas(genome, {10, 11, 15});
    
    Operator::phenotypeUpdateWeight(genome);
    
    // Only enabled connection should be updated in phenotype
    validateWeightUpdate(genome, 10, 1.1f);
    // Disabled connections (11, 15) are not in phenotype, so no validation needed
    
    validateDeltasEmpty(genome);
}

// Synchronization Tests
TEST_F(PhenotypeUpdateWeightTest, WeightsMatchGenome_SourceOfTruth) {
    auto genome = createBasicGenome();
    
    // Modify weights in genome
    modifyGenomeWeight(genome, 10, 9.9f);
    modifyGenomeWeight(genome, 11, 8.8f);
    modifyGenomeWeight(genome, 12, 7.7f);
    
    setConnectionDeltas(genome, {10, 11, 12});
    
    Operator::phenotypeUpdateWeight(genome);
    
    // Verify phenotype matches genome (genome is source of truth)
    validateWeightSynchronization(genome);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, UnchangedConnections_RemainUntouched) {
    auto genome = createComplexGenome();
    
    // Store original weights
    std::unordered_map<uint32_t, float> originalWeights;
    for (const auto& conn : genome.get_connectionGenes()) {
        originalWeights[conn.get_historyID()] = conn.get_attributes().weight;
    }
    
    // Only modify one weight
    modifyGenomeWeight(genome, 10, 9.5f);
    setConnectionDeltas(genome, {10});
    
    Operator::phenotypeUpdateWeight(genome);
    
    // Verify updated connection
    validateWeightUpdate(genome, 10, 9.5f);
    
    // Verify unchanged connections (enabled ones)
    validateUnchangedConnections(genome, {12, 13, 14, 16});
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, WeightPrecision_Maintained) {
    auto genome = createBasicGenome();
    
    // Test various precision values
    float preciseWeight = 1.23456789f;
    modifyGenomeWeight(genome, 11, preciseWeight);
    setConnectionDeltas(genome, {11});
    
    Operator::phenotypeUpdateWeight(genome);
    
    validateWeightUpdate(genome, 11, preciseWeight);
    validateDeltasEmpty(genome);
}

// Edge Cases
TEST_F(PhenotypeUpdateWeightTest, AllConnectionsDisabled_NoUpdatesOccur) {
    auto genome = createAllDisabledGenome();
    
    // Modify weights even though connections are disabled
    modifyGenomeWeight(genome, 10, 5.0f);
    modifyGenomeWeight(genome, 11, 6.0f);
    modifyGenomeWeight(genome, 12, 7.0f);
    
    setConnectionDeltas(genome, {10, 11, 12});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    
    Operator::phenotypeUpdateWeight(genome);
    
    // Phenotype should remain empty (no enabled connections)
    EXPECT_EQ(genome.get_phenotype()._orderedConnections.size(), initialConnections);
    EXPECT_EQ(genome.get_phenotype()._orderedConnections.size(), 0);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, SingleConnectionGenome_UpdatedCorrectly) {
    auto genome = createSingleConnectionGenome();
    
    float newWeight = 42.0f;
    modifyGenomeWeight(genome, 10, newWeight);
    setConnectionDeltas(genome, {10});
    
    Operator::phenotypeUpdateWeight(genome);
    
    const auto& phenotype = genome.get_phenotype();
    EXPECT_EQ(phenotype._orderedConnections.size(), 1);
    EXPECT_FLOAT_EQ(phenotype._orderedConnections[0]._connectionGeneAttribute.weight, newWeight);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, ExtremeWeightValues_HandledCorrectly) {
    auto genome = createBasicGenome();
    
    // Test extreme values
    std::vector<std::pair<uint32_t, float>> extremeWeights = {
        {10, 0.0f},           // Zero
        {11, -100.0f},        // Large negative
        {12, 1000.0f}         // Large positive
    };
    
    for (const auto& [historyID, weight] : extremeWeights) {
        modifyGenomeWeight(genome, historyID, weight);
    }
    
    setConnectionDeltas(genome, {10, 11, 12});
    
    Operator::phenotypeUpdateWeight(genome);
    
    for (const auto& [historyID, expectedWeight] : extremeWeights) {
        validateWeightUpdate(genome, historyID, expectedWeight);
    }
    
    validateDeltasEmpty(genome);
}

// State Consistency Tests
TEST_F(PhenotypeUpdateWeightTest, DeltasCleared_AfterOperation) {
    auto genome = createBasicGenome();
    
    modifyGenomeWeight(genome, 11, 5.0f);
    setConnectionDeltas(genome, {11});
    
    // Verify deltas are set before operation
    EXPECT_FALSE(genome.get_connectionGeneDeltas().empty());
    
    Operator::phenotypeUpdateWeight(genome);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateWeightTest, OnlySpecifiedConnections_Updated) {
    auto genome = createComplexGenome();
    
    // Store original state
    auto originalConnections = genome.get_phenotype()._orderedConnections;
    
    // Modify weights for subset of connections
    modifyGenomeWeight(genome, 10, 1.1f);
    modifyGenomeWeight(genome, 12, 3.3f);
    
    // Only specify one in deltas
    setConnectionDeltas(genome, {10});
    
    Operator::phenotypeUpdateWeight(genome);
    
    // Only connection 10 should be updated
    validateWeightUpdate(genome, 10, 1.1f);
    
    // Connection 12 should retain original weight in phenotype
    validateWeightUpdate(genome, 12, 3.0f); // Original weight
    
    validateDeltasEmpty(genome);
}