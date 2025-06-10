#include <gtest/gtest.h>
#include <unordered_map>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/operator/PhenotypeUpdateNode.hpp"

class PhenotypeUpdateNodeTest : public ::testing::Test {
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
        
        // Connections: B0->H2, I1->H2, H2->O3
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

    // Helper to create complex genome with multiple connections
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
        
        // Connections with direct I->O connection to split
        params._connectionHistoryIDs = {10, 11, 12, 13, 14, 15};
        params._sourceNodeHistoryIDs = {0, 1, 2, 3, 4, 1};
        params._targetNodeHistoryIDs = {3, 3, 4, 5, 6, 5};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},   // B0->H3
            ConnectionGeneAttributes{2.0f, true},   // I1->H3
            ConnectionGeneAttributes{3.0f, true},   // I2->H4
            ConnectionGeneAttributes{4.0f, true},   // H3->O5
            ConnectionGeneAttributes{5.0f, true},   // H4->O6
            ConnectionGeneAttributes{6.0f, true}    // I1->O5 (direct connection to split)
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to simulate node mutation: split connection by adding new node and connections
    void simulateNodeMutation(Genome& genome, uint32_t originalConnectionID, uint32_t newNodeID, uint32_t firstConnectionID, uint32_t secondConnectionID) {
        auto& connectionGenes = genome.get_connectionGenes();
        auto& nodeGenes = genome.get_nodeGenes();
        
        // Find and disable the original connection
        for (auto& conn : connectionGenes) {
            if (conn.get_historyID() == originalConnectionID) {
                conn.get_attributes().enabled = false;
                
                // Get source and target for the new connections
                const NodeGene& sourceNode = conn.get_sourceNodeGene();
                const NodeGene& targetNode = conn.get_targetNodeGene();
                
                // Add new hidden node
                nodeGenes.emplace_back(newNodeID, NodeType::HIDDEN, NodeGeneAttributes{ActivationType::SIGMOID});
                const NodeGene& newNode = nodeGenes.back();
                
                // Add first connection: source -> new node
                ConnectionGeneAttributes firstAttrs{1.0f, true};
                connectionGenes.emplace_back(firstConnectionID, sourceNode, newNode, firstAttrs);
                
                // Add second connection: new node -> target
                ConnectionGeneAttributes secondAttrs{conn.get_attributes().weight, true};
                connectionGenes.emplace_back(secondConnectionID, newNode, targetNode, secondAttrs);
                
                break;
            }
        }
    }

    // Helper to set connection deltas for node mutation pattern
    void setNodeMutationDeltas(Genome& genome, uint32_t disabledID, uint32_t firstNewID, uint32_t secondNewID) {
        auto& connectionDeltas = genome.get_connectionGeneDeltas();
        connectionDeltas = {disabledID, firstNewID, secondNewID};
    }

    // Validation helpers
    void validatePhenotypeStructure(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        
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
    }

    void validateNodeMutationPattern(const Genome& genome, uint32_t originalSourceID, uint32_t newNodeID, uint32_t originalTargetID) {
        const auto& phenotype = genome.get_phenotype();
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
        
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                includedHistoryIDs.insert(conn.get_sourceNodeGene().get_historyID());
                includedHistoryIDs.insert(conn.get_targetNodeGene().get_historyID());
            }
        }
        
        size_t phenotypeIndex = 0;
        for (const auto& node : nodeGenes) {
            if (includedHistoryIDs.count(node.get_historyID()) > 0) {
                historyToIndex[node.get_historyID()] = phenotypeIndex++;
            }
        }
        
        // Verify split pattern: originalSource -> newNode -> originalTarget
        bool foundFirstConnection = false;
        bool foundSecondConnection = false;
        
        for (const auto& phenConn : phenotype._orderedConnections) {
            // Check for originalSource -> newNode
            if (phenConn._sourceNodeIndex == historyToIndex[originalSourceID] &&
                phenConn._targetNodeIndex == historyToIndex[newNodeID]) {
                foundFirstConnection = true;
                EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
            }
            
            // Check for newNode -> originalTarget
            if (phenConn._sourceNodeIndex == historyToIndex[newNodeID] &&
                phenConn._targetNodeIndex == historyToIndex[originalTargetID]) {
                foundSecondConnection = true;
                EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
            }
        }
        
        EXPECT_TRUE(foundFirstConnection);
        EXPECT_TRUE(foundSecondConnection);
        
        // Verify new node is included in phenotype
        EXPECT_NE(historyToIndex.find(newNodeID), historyToIndex.end());
    }

    void validateCompleteReconstruction(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        const auto& nodeGenes = genome.get_nodeGenes();
        
        // Verify all enabled connections are represented
        size_t enabledCount = 0;
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                enabledCount++;
            }
        }
        
        EXPECT_EQ(phenotype._orderedConnections.size(), enabledCount);
        
        // Verify node inclusion logic
        std::unordered_set<uint32_t> expectedIncludedNodes;
        
        for (const auto& node : nodeGenes) {
            if (node.get_type() == NodeType::INPUT || 
                node.get_type() == NodeType::OUTPUT || 
                node.get_type() == NodeType::BIAS) {
                expectedIncludedNodes.insert(node.get_historyID());
            }
        }
        
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                expectedIncludedNodes.insert(conn.get_sourceNodeGene().get_historyID());
                expectedIncludedNodes.insert(conn.get_targetNodeGene().get_historyID());
            }
        }
        
        EXPECT_EQ(phenotype._nodeGeneAttributes.size(), expectedIncludedNodes.size());
    }
};

// Input Validation Tests
TEST_F(PhenotypeUpdateNodeTest, ExactlyThreeDeltas_Succeeds) {
    auto genome = createBasicGenome();
    
    // Simulate node mutation: split I1->H2 connection
    simulateNodeMutation(genome, 11, 4, 13, 14); // Split connection 11, add node 4, connections 13,14
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    EXPECT_NO_THROW(Operator::phenotypeUpdateNode(genome));
    
    validateNodeMutationPattern(genome, 1, 4, 2);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateNodeTest, WrongNumberOfDeltas_Fails) {
    auto genome = createBasicGenome();
    
    // Test empty deltas
    setNodeMutationDeltas(genome, 0, 0, 0);
    genome.get_connectionGeneDeltas().clear();
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "Connection deltas must contain exactly 3 IDs");
    
    // Test too few deltas
    genome.get_connectionGeneDeltas() = {11};
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "Connection deltas must contain exactly 3 IDs");
    
    // Test too many deltas
    genome.get_connectionGeneDeltas() = {11, 13, 14, 15};
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "Connection deltas must contain exactly 3 IDs");
}

TEST_F(PhenotypeUpdateNodeTest, InvalidConnectionIds_Fails) {
    auto genome = createBasicGenome();
    
    // Test non-existent connection IDs
    setNodeMutationDeltas(genome, 999, 998, 997);
    
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "Disabled connection not found in genome");
}

TEST_F(PhenotypeUpdateNodeTest, IncorrectConnectionStates_Fails) {
    auto genome = createBasicGenome();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    
    // Test disabled connection that should be disabled - this should work
    setNodeMutationDeltas(genome, 11, 13, 14);
    EXPECT_NO_THROW(Operator::phenotypeUpdateNode(genome));
    
    // Reset for next test
    genome = createBasicGenome();
    simulateNodeMutation(genome, 11, 4, 13, 14);
    
    // Manually enable the "disabled" connection - should fail
    for (auto& conn : genome.get_connectionGenes()) {
        if (conn.get_historyID() == 11) {
            conn.get_attributes().enabled = true;
            break;
        }
    }
    setNodeMutationDeltas(genome, 11, 13, 14);
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "Disabled connection must be disabled");
    
    // Reset for next test
    genome = createBasicGenome();
    simulateNodeMutation(genome, 11, 4, 13, 14);
    
    // Manually disable one of the "enabled" connections - should fail
    for (auto& conn : genome.get_connectionGenes()) {
        if (conn.get_historyID() == 13) {
            conn.get_attributes().enabled = false;
            break;
        }
    }
    setNodeMutationDeltas(genome, 11, 13, 14);
    EXPECT_DEATH(Operator::phenotypeUpdateNode(genome), "First new connection must be enabled");
}

// Node Mutation Pattern Tests
TEST_F(PhenotypeUpdateNodeTest, FirstDelta_IsDisabledConnection) {
    auto genome = createBasicGenome();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    // Find the disabled connection
    const ConnectionGene* disabledConnection = nullptr;
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_historyID() == 11) {
            disabledConnection = &conn;
            break;
        }
    }
    
    ASSERT_NE(disabledConnection, nullptr);
    EXPECT_FALSE(disabledConnection->get_attributes().enabled);
    
    Operator::phenotypeUpdateNode(genome);
    
    // Verify the pattern was correctly processed
    validateNodeMutationPattern(genome, 1, 4, 2);
}

TEST_F(PhenotypeUpdateNodeTest, SecondThirdDeltas_AreEnabledConnections) {
    auto genome = createBasicGenome();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    // Find the new connections
    const ConnectionGene* firstNewConnection = nullptr;
    const ConnectionGene* secondNewConnection = nullptr;
    
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_historyID() == 13) {
            firstNewConnection = &conn;
        } else if (conn.get_historyID() == 14) {
            secondNewConnection = &conn;
        }
    }
    
    ASSERT_NE(firstNewConnection, nullptr);
    ASSERT_NE(secondNewConnection, nullptr);
    EXPECT_TRUE(firstNewConnection->get_attributes().enabled);
    EXPECT_TRUE(secondNewConnection->get_attributes().enabled);
    
    Operator::phenotypeUpdateNode(genome);
    
    validateNodeMutationPattern(genome, 1, 4, 2);
}

TEST_F(PhenotypeUpdateNodeTest, SplitPattern_CorrectlyFormed) {
    auto genome = createBasicGenome();
    
    // Split I1->H2 connection into I1->H4->H2
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    Operator::phenotypeUpdateNode(genome);
    
    validateNodeMutationPattern(genome, 1, 4, 2);
    validateCompleteReconstruction(genome);
}

TEST_F(PhenotypeUpdateNodeTest, SplitNode_IncludedInPhenotype) {
    auto genome = createBasicGenome();
    
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should have one more node
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 1);
    
    validateNodeMutationPattern(genome, 1, 4, 2);
}

// Phenotype Reconstruction Tests
TEST_F(PhenotypeUpdateNodeTest, AllNodes_ReEvaluatedForInclusion) {
    auto genome = createComplexGenome();
    
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    // Split a connection
    simulateNodeMutation(genome, 15, 7, 16, 17); // Split I1->O5
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should have one more node (the new hidden node)
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 1);
    
    validateCompleteReconstruction(genome);
    validateNodeMutationPattern(genome, 1, 7, 5);
}

TEST_F(PhenotypeUpdateNodeTest, NewHiddenNode_IncludedWithConnections) {
    auto genome = createBasicGenome();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& phenotype = genome.get_phenotype();
    const auto& nodeGenes = genome.get_nodeGenes();
    
    // Find new node in phenotype
    bool foundNewNode = false;
    for (const auto& node : nodeGenes) {
        if (node.get_historyID() == 4 && node.get_type() == NodeType::HIDDEN) {
            foundNewNode = true;
            break;
        }
    }
    EXPECT_TRUE(foundNewNode);
    
    validateNodeMutationPattern(genome, 1, 4, 2);
}

TEST_F(PhenotypeUpdateNodeTest, NetworkConnectivity_ProperlyUpdated) {
    auto genome = createComplexGenome();
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    
    // Split a connection - this disables 1 and adds 2, so net +1
    simulateNodeMutation(genome, 15, 7, 16, 17);
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should have one more connection (disabled 1, added 2)
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections + 1);
    
    validateCompleteReconstruction(genome);
}

TEST_F(PhenotypeUpdateNodeTest, AllReferences_RemainValid) {
    auto genome = createComplexGenome();
    
    simulateNodeMutation(genome, 15, 7, 16, 17);
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    validatePhenotypeStructure(genome);
    validateCompleteReconstruction(genome);
}

// Edge Cases
TEST_F(PhenotypeUpdateNodeTest, SplitCreatingFirstHiddenNode) {
    // Create genome with only I/O/B nodes
    GenomeParams params;
    
    params._nodeHistoryIDs = {0, 1, 2};
    params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::OUTPUT};
    params._nodeAttributes = {
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID},
        NodeGeneAttributes{ActivationType::SIGMOID}
    };
    
    params._connectionHistoryIDs = {10, 11};
    params._sourceNodeHistoryIDs = {0, 1};
    params._targetNodeHistoryIDs = {2, 2};
    params._connectionAttributes = {
        ConnectionGeneAttributes{1.0f, true},
        ConnectionGeneAttributes{2.0f, true}
    };
    
    Genome genome(params);
    Operator::phenotypeConstruct(genome);
    
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    // Split I1->O2 to create first hidden node
    simulateNodeMutation(genome, 11, 3, 12, 13);
    setNodeMutationDeltas(genome, 11, 12, 13);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should have one more node (first hidden)
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 1);
    
    validateNodeMutationPattern(genome, 1, 3, 2);
}

TEST_F(PhenotypeUpdateNodeTest, SplitInputToOutputConnection) {
    auto genome = createComplexGenome();
    
    // Split direct I1->O5 connection
    simulateNodeMutation(genome, 15, 7, 16, 17);
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    validateNodeMutationPattern(genome, 1, 7, 5);
    validateCompleteReconstruction(genome);
}

TEST_F(PhenotypeUpdateNodeTest, SplitInComplexNetwork_MaintainsOtherConnections) {
    auto genome = createComplexGenome();
    
    // Count enabled connections before
    size_t enabledBefore = 0;
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_attributes().enabled) {
            enabledBefore++;
        }
    }
    
    simulateNodeMutation(genome, 15, 7, 16, 17);
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    // Should have disabled 1, added 2, so net +1
    size_t enabledAfter = 0;
    for (const auto& conn : genome.get_connectionGenes()) {
        if (conn.get_attributes().enabled) {
            enabledAfter++;
        }
    }
    
    EXPECT_EQ(enabledAfter, enabledBefore + 1);
    
    validateCompleteReconstruction(genome);
}

// State Consistency Tests
TEST_F(PhenotypeUpdateNodeTest, DeltasCleared_AfterOperation) {
    auto genome = createBasicGenome();
    
    simulateNodeMutation(genome, 11, 4, 13, 14);
    setNodeMutationDeltas(genome, 11, 13, 14);
    
    Operator::phenotypeUpdateNode(genome);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateNodeTest, PhenotypeCompletelyRebuilt) {
    auto genome = createComplexGenome();
    
    // Store initial state
    auto initialPhenotype = genome.get_phenotype();
    
    simulateNodeMutation(genome, 15, 7, 16, 17);
    setNodeMutationDeltas(genome, 15, 16, 17);
    
    Operator::phenotypeUpdateNode(genome);
    
    const auto& newPhenotype = genome.get_phenotype();
    
    // Phenotype should be different (completely rebuilt)
    EXPECT_NE(initialPhenotype._nodeGeneAttributes.size(), newPhenotype._nodeGeneAttributes.size());
    
    validateCompleteReconstruction(genome);
    validatePhenotypeStructure(genome);
}