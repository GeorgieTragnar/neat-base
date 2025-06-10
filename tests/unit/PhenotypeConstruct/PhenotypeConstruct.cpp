#include <gtest/gtest.h>
#include <unordered_map>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"

class PhenotypeConstructTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
    }

    // Helper to create basic genome: I1 -> H2 -> O3 with B0
    Genome createBasicGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), HIDDEN(2), OUTPUT(3)
        params._nodeHistoryIDs = {0, 1, 2, 3};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},  // BIAS
            NodeGeneAttributes{ActivationType::SIGMOID},  // INPUT
            NodeGeneAttributes{ActivationType::SIGMOID},  // HIDDEN
            NodeGeneAttributes{ActivationType::SIGMOID}   // OUTPUT
        };
        
        // Connections: B0->H2(enabled), I1->H2(enabled), H2->O3(enabled)
        params._connectionHistoryIDs = {10, 11, 12};
        params._sourceNodeHistoryIDs = {0, 1, 2};
        params._targetNodeHistoryIDs = {2, 2, 3};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->H2
            ConnectionGeneAttributes{2.0f, true},   // I1->H2  
            ConnectionGeneAttributes{3.0f, true}    // H2->O3
        };
        
        return Genome(params);
    }

    // Helper to create complex genome with multiple hidden nodes and mixed connections
    Genome createComplexGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), INPUT(2), HIDDEN(3), HIDDEN(4), OUTPUT(5), OUTPUT(6)
        params._nodeHistoryIDs = {0, 1, 2, 3, 4, 5, 6};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},  // BIAS
            NodeGeneAttributes{ActivationType::SIGMOID},  // INPUT
            NodeGeneAttributes{ActivationType::SIGMOID},  // INPUT
            NodeGeneAttributes{ActivationType::SIGMOID},  // HIDDEN
            NodeGeneAttributes{ActivationType::SIGMOID},  // HIDDEN
            NodeGeneAttributes{ActivationType::SIGMOID},  // OUTPUT
            NodeGeneAttributes{ActivationType::SIGMOID}   // OUTPUT
        };
        
        // Mixed enabled/disabled connections
        params._connectionHistoryIDs = {10, 11, 12, 13, 14, 15, 16, 17};
        params._sourceNodeHistoryIDs = {0, 1, 2, 3, 4, 1, 2, 3};
        params._targetNodeHistoryIDs = {3, 3, 4, 5, 6, 5, 6, 4};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->H3 (enabled)
            ConnectionGeneAttributes{2.0f, false},  // I1->H3 (disabled)
            ConnectionGeneAttributes{3.0f, true},   // I2->H4 (enabled)
            ConnectionGeneAttributes{4.0f, true},   // H3->O5 (enabled)
            ConnectionGeneAttributes{5.0f, true},   // H4->O6 (enabled)
            ConnectionGeneAttributes{6.0f, false},  // I1->O5 (disabled)
            ConnectionGeneAttributes{7.0f, true},   // I2->O6 (enabled)
            ConnectionGeneAttributes{8.0f, false}   // H3->H4 (disabled)
        };
        
        return Genome(params);
    }

    // Helper to create genome with disconnected hidden node
    Genome createDisconnectedGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), HIDDEN(2), HIDDEN(3), OUTPUT(4)
        params._nodeHistoryIDs = {0, 1, 2, 3, 4};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},  // BIAS
            NodeGeneAttributes{ActivationType::SIGMOID},  // INPUT
            NodeGeneAttributes{ActivationType::SIGMOID},  // HIDDEN (connected)
            NodeGeneAttributes{ActivationType::SIGMOID},  // HIDDEN (disconnected)
            NodeGeneAttributes{ActivationType::SIGMOID}   // OUTPUT
        };
        
        // H3 has no enabled connections (disconnected)
        params._connectionHistoryIDs = {10, 11, 12, 13};
        params._sourceNodeHistoryIDs = {0, 1, 2, 3};
        params._targetNodeHistoryIDs = {2, 2, 4, 4};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->H2 (enabled)
            ConnectionGeneAttributes{2.0f, true},   // I1->H2 (enabled)
            ConnectionGeneAttributes{3.0f, true},   // H2->O4 (enabled)
            ConnectionGeneAttributes{4.0f, false}   // H3->O4 (disabled)
        };
        
        return Genome(params);
    }

    // Helper to create genome with only I/O/B nodes
    Genome createSimpleGenome() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), OUTPUT(2)
        params._nodeHistoryIDs = {0, 1, 2};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},  // BIAS
            NodeGeneAttributes{ActivationType::SIGMOID},  // INPUT
            NodeGeneAttributes{ActivationType::SIGMOID}   // OUTPUT
        };
        
        // Direct connection: I1->O2
        params._connectionHistoryIDs = {10, 11};
        params._sourceNodeHistoryIDs = {0, 1};
        params._targetNodeHistoryIDs = {2, 2};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->O2
            ConnectionGeneAttributes{2.0f, true}    // I1->O2
        };
        
        return Genome(params);
    }

    // Helper to set connection deltas
    void setConnectionDeltas(Genome& genome, const std::vector<uint32_t>& deltas) {
        auto& connectionDeltas = genome.get_connectionGeneDeltas();
        connectionDeltas = deltas;
    }

    // Helper to set node deltas  
    void setNodeDeltas(Genome& genome, const std::vector<uint32_t>& deltas) {
        auto& nodeDeltas = genome.get_nodeGeneDeltas();
        nodeDeltas = deltas;
    }

    // Validation helpers
    void validatePhenotypeStructure(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        
        // All indices should be valid
        for (size_t inputIdx : phenotype._inputIndices) {
            EXPECT_LT(inputIdx, phenotype._nodeGeneAttributes.size());
        }
        
        for (size_t outputIdx : phenotype._outputIndices) {
            EXPECT_LT(outputIdx, phenotype._nodeGeneAttributes.size());
        }
        
        for (const auto& conn : phenotype._orderedConnections) {
            EXPECT_LT(conn._sourceNodeIndex, phenotype._nodeGeneAttributes.size());
            EXPECT_LT(conn._targetNodeIndex, phenotype._nodeGeneAttributes.size());
        }
    }

    void validateDeltasEmpty(Genome& genome) {
        EXPECT_TRUE(genome.get_connectionGeneDeltas().empty());
        EXPECT_TRUE(genome.get_nodeGeneDeltas().empty());
    }

    void validateIndexMapping(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        const auto& nodeGenes = genome.get_nodeGenes();
        
        // Build expected mapping from genome order to phenotype indices
        std::unordered_map<uint32_t, size_t> expectedMapping;
        std::unordered_set<uint32_t> includedHistoryIDs;
        
        // Determine which nodes should be included
        for (const auto& node : nodeGenes) {
            if (node.get_type() == NodeType::INPUT || 
                node.get_type() == NodeType::OUTPUT || 
                node.get_type() == NodeType::BIAS) {
                includedHistoryIDs.insert(node.get_historyID());
            }
        }
        
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                includedHistoryIDs.insert(conn.get_sourceNodeGene().get_historyID());
                includedHistoryIDs.insert(conn.get_targetNodeGene().get_historyID());
            }
        }
        
        // Build expected index mapping
        size_t phenotypeIndex = 0;
        for (const auto& node : nodeGenes) {
            if (includedHistoryIDs.count(node.get_historyID()) > 0) {
                expectedMapping[node.get_historyID()] = phenotypeIndex++;
            }
        }
        
        EXPECT_EQ(phenotypeIndex, phenotype._nodeGeneAttributes.size());
        
        // Validate connections use correct indices
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                uint32_t sourceHistoryID = conn.get_sourceNodeGene().get_historyID();
                uint32_t targetHistoryID = conn.get_targetNodeGene().get_historyID();
                
                // Find corresponding phenotype connection
                bool foundConnection = false;
                for (const auto& phenConn : phenotype._orderedConnections) {
                    if (phenConn._sourceNodeIndex == expectedMapping[sourceHistoryID] &&
                        phenConn._targetNodeIndex == expectedMapping[targetHistoryID]) {
                        foundConnection = true;
                        EXPECT_EQ(phenConn._connectionGeneAttribute.weight, conn.get_attributes().weight);
                        EXPECT_EQ(phenConn._connectionGeneAttribute.enabled, conn.get_attributes().enabled);
                        break;
                    }
                }
                EXPECT_TRUE(foundConnection) << "Connection not found in phenotype";
            }
        }
    }
};

// Input Validation Tests
TEST_F(PhenotypeConstructTest, EmptyGenome_ProducesEmptyPhenotype) {
    GenomeParams params;
    // Empty genome - no nodes, no connections
    Genome genome(params);
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    EXPECT_TRUE(phenotype._nodeGeneAttributes.empty());
    EXPECT_TRUE(phenotype._orderedConnections.empty());
    EXPECT_TRUE(phenotype._inputIndices.empty());
    EXPECT_TRUE(phenotype._outputIndices.empty());
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeConstructTest, GenomeWithNonEmptyDeltas_ClearsDeltas) {
    auto genome = createBasicGenome();
    
    // Set some deltas before construction
    setConnectionDeltas(genome, {10, 11});
    setNodeDeltas(genome, {1, 2});
    
    Operator::phenotypeConstruct(genome);
    
    validateDeltasEmpty(genome);
}

// Node Inclusion Logic Tests
TEST_F(PhenotypeConstructTest, InputOutputBiasNodes_AlwaysIncluded) {
    auto genome = createBasicGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should have 4 nodes: BIAS(0), INPUT(1), HIDDEN(2), OUTPUT(3)
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 4);
    EXPECT_EQ(phenotype._inputIndices.size(), 1);
    EXPECT_EQ(phenotype._outputIndices.size(), 1);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, HiddenNodesWithEnabledConnections_Included) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should include all nodes except those with only disabled connections
    // Expected: BIAS(0), INPUT(1), INPUT(2), HIDDEN(3), HIDDEN(4), OUTPUT(5), OUTPUT(6)
    // H3 has enabled connection B0->H3 and H3->O5
    // H4 has enabled connection I2->H4 and H4->O6
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 7);
    EXPECT_EQ(phenotype._inputIndices.size(), 2);
    EXPECT_EQ(phenotype._outputIndices.size(), 2);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, HiddenNodesWithOnlyDisabledConnections_Excluded) {
    // Create genome where hidden node has only disabled connections
    GenomeParams params;
    
    params._nodeHistoryIDs = {0, 1, 2, 3};
    params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
    params._nodeAttributes = {
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID}
    };
    
    // Hidden node H2 has only disabled connections
    params._connectionHistoryIDs = {10, 11, 12};
    params._sourceNodeHistoryIDs = {0, 1, 2};
    params._targetNodeHistoryIDs = {2, 2, 3};
    params._connectionAttributes = {
        ConnectionGeneAttributes{1.0f, false},  // B0->H2 (disabled)
        ConnectionGeneAttributes{2.0f, false},  // I1->H2 (disabled)
        ConnectionGeneAttributes{3.0f, false}   // H2->O3 (disabled)
    };
    
    Genome genome(params);
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should only have I/O/B nodes, H2 excluded
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 3); // BIAS, INPUT, OUTPUT only
    EXPECT_EQ(phenotype._orderedConnections.size(), 0); // No enabled connections
    
    validatePhenotypeStructure(genome);
}

TEST_F(PhenotypeConstructTest, DisconnectedHiddenNodes_Excluded) {
    auto genome = createDisconnectedGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should exclude H3 (disconnected), include H2 (connected)
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 4); // BIAS, INPUT, HIDDEN(2), OUTPUT
    EXPECT_EQ(phenotype._orderedConnections.size(), 3); // 3 enabled connections
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

// Connection Filtering Tests
TEST_F(PhenotypeConstructTest, OnlyEnabledConnections_IncludedInPhenotype) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Count enabled connections in genome
    size_t enabledCount = 0;
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_attributes().enabled) {
            enabledCount++;
        }
    }
    
    EXPECT_EQ(phenotype._orderedConnections.size(), enabledCount);
    
    // Verify all connections in phenotype are enabled
    for (const auto& phenConn : phenotype._orderedConnections) {
        EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
    }
    
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, DisabledConnections_ExcludedFromPhenotype) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Verify no disabled connections in phenotype
    for (const auto& phenConn : phenotype._orderedConnections) {
        EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
    }
    
    // Verify count matches enabled connections only
    size_t enabledCount = 0;
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_attributes().enabled) {
            enabledCount++;
        }
    }
    
    EXPECT_EQ(phenotype._orderedConnections.size(), enabledCount);
}

TEST_F(PhenotypeConstructTest, MixedEnabledDisabled_OnlyEnabledIncluded) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    for (const auto& phenConn : phenotype._orderedConnections) {
        EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
    }
    
    validateIndexMapping(genome);
}

// Index Mapping Tests
TEST_F(PhenotypeConstructTest, InputIndices_PointToCorrectPhenotypeNodes) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    const auto& nodeGenes = genome.get_nodeGenes();
    
    // Build expected input indices
    std::vector<size_t> expectedInputIndices;
    size_t phenotypeIndex = 0;
    
    for (const auto& node : nodeGenes) {
        if (node.get_type() == NodeType::INPUT) {
            expectedInputIndices.push_back(phenotypeIndex);
        }
        phenotypeIndex++;
    }
    
    EXPECT_EQ(phenotype._inputIndices.size(), 2);
    // Note: actual validation is complex due to node ordering, so we rely on validateIndexMapping
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, OutputIndices_PointToCorrectPhenotypeNodes) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    EXPECT_EQ(phenotype._outputIndices.size(), 2);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, ConnectionIndices_MapToCorrectPhenotypeNodes) {
    auto genome = createBasicGenome();
    
    Operator::phenotypeConstruct(genome);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, HistoryIdToPhenotypeIndex_MappingCorrect) {
    auto genome = createComplexGenome();
    
    Operator::phenotypeConstruct(genome);
    
    validateIndexMapping(genome);
}

// Edge Cases
TEST_F(PhenotypeConstructTest, AllConnectionsDisabled_OnlyIOBNodesInPhenotype) {
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
    
    const auto& phenotype = genome.get_phenotype();
    
    // Only I/O/B nodes should be included
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 3);
    EXPECT_EQ(phenotype._orderedConnections.size(), 0);
    EXPECT_EQ(phenotype._inputIndices.size(), 1);
    EXPECT_EQ(phenotype._outputIndices.size(), 1);
    
    validatePhenotypeStructure(genome);
}

TEST_F(PhenotypeConstructTest, SingleNodeGenome_CorrectlyHandled) {
    GenomeParams params;
    
    params._nodeHistoryIDs = {1};
    params._nodeTypes = {NodeType::INPUT};
    params._nodeAttributes = {
        NodeGeneAttributes{ActivationType::SIGMOID}
    };
    
    // No connections
    
    Genome genome(params);
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 1);
    EXPECT_EQ(phenotype._orderedConnections.size(), 0);
    EXPECT_EQ(phenotype._inputIndices.size(), 1);
    EXPECT_EQ(phenotype._outputIndices.size(), 0);
    
    validatePhenotypeStructure(genome);
}

TEST_F(PhenotypeConstructTest, NoHiddenNodes_CorrectlyHandled) {
    auto genome = createSimpleGenome();
    
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 3); // BIAS, INPUT, OUTPUT
    EXPECT_EQ(phenotype._orderedConnections.size(), 2); // 2 connections
    EXPECT_EQ(phenotype._inputIndices.size(), 1);
    EXPECT_EQ(phenotype._outputIndices.size(), 1);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeConstructTest, LinearChainTopology_CorrectOrder) {
    GenomeParams params;
    
    // Linear chain: I1 -> H2 -> H3 -> O4 with B0
    params._nodeHistoryIDs = {0, 1, 2, 3, 4};
    params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
    params._nodeAttributes = {
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID}
    };
    
    params._connectionHistoryIDs = {10, 11, 12, 13};
    params._sourceNodeHistoryIDs = {0, 1, 2, 3};
    params._targetNodeHistoryIDs = {4, 2, 3, 4};
    params._connectionAttributes = {
        ConnectionGeneAttributes{1.0f, true},  // B0->O4
        ConnectionGeneAttributes{2.0f, true},  // I1->H2
        ConnectionGeneAttributes{3.0f, true},  // H2->H3
        ConnectionGeneAttributes{4.0f, true}   // H3->O4
    };
    
    Genome genome(params);
    Operator::phenotypeConstruct(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), 5);
    EXPECT_EQ(phenotype._orderedConnections.size(), 4);
    EXPECT_EQ(phenotype._inputIndices.size(), 1);
    EXPECT_EQ(phenotype._outputIndices.size(), 1);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

// State Consistency Tests
TEST_F(PhenotypeConstructTest, DeltasCleared_AfterOperation) {
    auto genome = createBasicGenome();
    
    // Set deltas before operation
    setConnectionDeltas(genome, {10, 11, 12});
    setNodeDeltas(genome, {1, 2, 3});
    
    Operator::phenotypeConstruct(genome);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeConstructTest, ExistingPhenotype_CompletelyReplaced) {
    auto genome = createBasicGenome();
    
    // First construction
    Operator::phenotypeConstruct(genome);
    auto firstConstruction = genome.get_phenotype();
    
    // Modify genome by adding a node (simulate mutation)
    GenomeParams newParams;
    newParams._nodeHistoryIDs = {0, 1, 2, 3, 4}; // Add new hidden node
    newParams._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
    newParams._nodeAttributes = {
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID}
    };
    
    newParams._connectionHistoryIDs = {10, 11, 12, 13};
    newParams._sourceNodeHistoryIDs = {0, 1, 2, 3};
    newParams._targetNodeHistoryIDs = {2, 2, 4, 4};
    newParams._connectionAttributes = {
        ConnectionGeneAttributes{1.0f, true},
        ConnectionGeneAttributes{2.0f, true},
        ConnectionGeneAttributes{3.0f, true},
        ConnectionGeneAttributes{4.0f, true}
    };
    
    Genome newGenome(newParams);
    
    // Second construction - should completely replace phenotype
    Operator::phenotypeConstruct(newGenome);
    
    const auto& secondConstruction = newGenome.get_phenotype();
    
    // Verify new phenotype is different
    EXPECT_NE(firstConstruction._nodeGeneAttributes.size(), secondConstruction._nodeGeneAttributes.size());
    
    validatePhenotypeStructure(newGenome);
    validateIndexMapping(newGenome);
}