#include <gtest/gtest.h>
#include <unordered_map>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/operator/PhenotypeUpdateConnection.hpp"

class PhenotypeUpdateConnectionTest : public ::testing::Test {
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
        Operator::phenotypeConstruct(genome); // Build initial phenotype
        return genome;
    }

    // Helper to create genome with extra hidden node but no connections to it
    Genome createGenomeWithUnconnectedHidden() {
        GenomeParams params;
        
        // Nodes: BIAS(0), INPUT(1), HIDDEN(2), HIDDEN(4), OUTPUT(3)
        params._nodeHistoryIDs = {0, 1, 2, 4, 3};
        params._nodeTypes = {NodeType::BIAS, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID},
            NodeGeneAttributes{ActivationType::SIGMOID}
        };
        
        // Connections: B0->H2, I1->H2, H2->O3 (H4 not connected)
        params._connectionHistoryIDs = {10, 11, 12};
        params._sourceNodeHistoryIDs = {0, 1, 2};
        params._targetNodeHistoryIDs = {2, 2, 3};
        params._connectionAttributes = {
            ConnectionGeneAttributes{1.0f, true},
            ConnectionGeneAttributes{2.0f, true},
            ConnectionGeneAttributes{3.0f, true}
        };
        
        Genome genome(params);
        Operator::phenotypeConstruct(genome); // Build initial phenotype (won't include H4)
        return genome;
    }

    // Helper to create empty genome
    Genome createEmptyGenome() {
        GenomeParams params;
        Genome genome(params);
        Operator::phenotypeConstruct(genome);
        return genome;
    }

    // Helper to add connection to genome after initial construction
    void addConnectionToGenome(Genome& genome, uint32_t historyID, uint32_t sourceHistoryID, uint32_t targetHistoryID, float weight = 5.0f, bool enabled = true) {
        // This is a simplified helper - in real code, you'd use proper mutation operators
        // For testing, we'll manually add to the genome structure
        auto& connectionGenes = genome.get_connectionGenes();
        auto& nodeGenes = genome.get_nodeGenes();
        
        // Find source and target nodes
        const NodeGene* sourceNode = nullptr;
        const NodeGene* targetNode = nullptr;
        
        for (const auto& node : nodeGenes) {
            if (node.get_historyID() == sourceHistoryID) {
                sourceNode = &node;
            }
            if (node.get_historyID() == targetHistoryID) {
                targetNode = &node;
            }
        }
        
        if (sourceNode && targetNode) {
            ConnectionGeneAttributes attrs{weight, enabled};
            ConnectionGene newConnection(historyID, *sourceNode, *targetNode, attrs);
            connectionGenes.push_back(newConnection);
        }
    }

    // Helper to set connection deltas
    void setConnectionDeltas(Genome& genome, const std::vector<uint32_t>& deltas) {
        auto& connectionDeltas = genome.get_connectionGeneDeltas();
        connectionDeltas = deltas;
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

    void validateIndexMapping(const Genome& genome) {
        const auto& phenotype = genome.get_phenotype();
        const auto& nodeGenes = genome.get_nodeGenes();
        
        // Build expected mapping
        std::unordered_map<uint32_t, size_t> expectedMapping;
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
                expectedMapping[node.get_historyID()] = phenotypeIndex++;
            }
        }
        
        EXPECT_EQ(phenotypeIndex, phenotype._nodeGeneAttributes.size());
        
        // Validate connections use correct indices
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                uint32_t sourceHistoryID = conn.get_sourceNodeGene().get_historyID();
                uint32_t targetHistoryID = conn.get_targetNodeGene().get_historyID();
                
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
                EXPECT_TRUE(foundConnection);
            }
        }
    }

    void validateNewConnectionExists(const Genome& genome, uint32_t sourceHistoryID, uint32_t targetHistoryID, float expectedWeight) {
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
        
        // Find the connection
        bool foundConnection = false;
        for (const auto& phenConn : phenotype._orderedConnections) {
            if (phenConn._sourceNodeIndex == historyToIndex[sourceHistoryID] &&
                phenConn._targetNodeIndex == historyToIndex[targetHistoryID]) {
                foundConnection = true;
                EXPECT_FLOAT_EQ(phenConn._connectionGeneAttribute.weight, expectedWeight);
                EXPECT_TRUE(phenConn._connectionGeneAttribute.enabled);
                break;
            }
        }
        EXPECT_TRUE(foundConnection);
    }
};

// Input Validation Tests
TEST_F(PhenotypeUpdateConnectionTest, ExactlyOneDelta_Succeeds) {
    auto genome = createBasicGenome();
    
    // Add a new connection to genome
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true); // I1->O3
    setConnectionDeltas(genome, {13});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    
    EXPECT_NO_THROW(Operator::phenotypeUpdateConnection(genome));
    
    // Should have one more connection
    EXPECT_EQ(genome.get_phenotype()._orderedConnections.size(), initialConnections + 1);
    validateNewConnectionExists(genome, 1, 3, 4.0f);
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, EmptyDeltas_Fails) {
    auto genome = createBasicGenome();
    
    setConnectionDeltas(genome, {});
    
    EXPECT_DEATH(Operator::phenotypeUpdateConnection(genome), "Connection deltas must contain exactly 1 ID");
}

TEST_F(PhenotypeUpdateConnectionTest, MultipleDelta_Fails) {
    auto genome = createBasicGenome();
    
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    addConnectionToGenome(genome, 14, 0, 3, 5.0f, true);
    setConnectionDeltas(genome, {13, 14});
    
    EXPECT_DEATH(Operator::phenotypeUpdateConnection(genome), "Connection deltas must contain exactly 1 ID");
}

TEST_F(PhenotypeUpdateConnectionTest, NonExistentConnectionId_Fails) {
    auto genome = createBasicGenome();
    
    setConnectionDeltas(genome, {999}); // Non-existent ID
    
    EXPECT_DEATH(Operator::phenotypeUpdateConnection(genome), "New connection with history ID not found in genome");
}

TEST_F(PhenotypeUpdateConnectionTest, DisabledConnection_Fails) {
    auto genome = createBasicGenome();
    
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, false); // Disabled connection
    setConnectionDeltas(genome, {13});
    
    EXPECT_DEATH(Operator::phenotypeUpdateConnection(genome), "New connection must be enabled");
}

// Connection Addition Tests
TEST_F(PhenotypeUpdateConnectionTest, ConnectionBetweenExistingNodes_AddedToPhenotype) {
    auto genome = createBasicGenome();
    
    // Add connection between existing nodes I1->O3
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should add connection but not change node count
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections + 1);
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes);
    
    validateNewConnectionExists(genome, 1, 3, 4.0f);
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, ConnectionRequiringSourceNode_AddsNodeAndConnection) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    // Add connection from unconnected H4 to existing O3
    addConnectionToGenome(genome, 13, 4, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should add both node and connection
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections + 1);
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 1);
    
    validateNewConnectionExists(genome, 4, 3, 4.0f);
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, ConnectionRequiringTargetNode_AddsNodeAndConnection) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    // Add connection from existing I1 to unconnected H4
    addConnectionToGenome(genome, 13, 1, 4, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should add both node and connection
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections + 1);
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 1);
    
    validateNewConnectionExists(genome, 1, 4, 4.0f);
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, ConnectionRequiringBothNodes_AddsBothNodesAndConnection) {
    auto genome = createBasicGenome();
    
    // Add two new hidden nodes
    auto& nodeGenes = genome.get_nodeGenes();
    nodeGenes.emplace_back(5, NodeType::HIDDEN, NodeGeneAttributes{ActivationType::SIGMOID});
    nodeGenes.emplace_back(6, NodeType::HIDDEN, NodeGeneAttributes{ActivationType::SIGMOID});
    
    // Add connection between the two new nodes
    addConnectionToGenome(genome, 13, 5, 6, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    size_t initialConnections = genome.get_phenotype()._orderedConnections.size();
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Should add both nodes and connection
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections + 1);
    EXPECT_EQ(phenotype._nodeGeneAttributes.size(), initialNodes + 2);
    
    validateNewConnectionExists(genome, 5, 6, 4.0f);
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

// Node Addition Logic Tests
TEST_F(PhenotypeUpdateConnectionTest, NewNodesAdded_CorrectAttributes) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    // Add connection to unconnected H4
    addConnectionToGenome(genome, 13, 4, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    const auto& nodeGenes = genome.get_nodeGenes();
    
    // Find H4 in genome
    const NodeGene* h4Node = nullptr;
    for (const auto& node : nodeGenes) {
        if (node.get_historyID() == 4) {
            h4Node = &node;
            break;
        }
    }
    ASSERT_NE(h4Node, nullptr);
    
    // Find H4 in phenotype
    bool foundH4InPhenotype = false;
    for (const auto& phenAttr : phenotype._nodeGeneAttributes) {
        if (phenAttr.activationType == h4Node->get_attributes().activationType) {
            foundH4InPhenotype = true;
            break;
        }
    }
    EXPECT_TRUE(foundH4InPhenotype);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, IndexMapping_UpdatedForNewNodes) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    addConnectionToGenome(genome, 13, 4, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, InputOutputIndices_RecalculatedWhenNeeded) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    size_t initialInputs = genome.get_phenotype()._inputIndices.size();
    size_t initialOutputs = genome.get_phenotype()._outputIndices.size();
    
    addConnectionToGenome(genome, 13, 4, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // Input/output counts should remain same
    EXPECT_EQ(phenotype._inputIndices.size(), initialInputs);
    EXPECT_EQ(phenotype._outputIndices.size(), initialOutputs);
    
    // But indices might have changed due to node reordering
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, ExistingConnections_RemainValid) {
    auto genome = createBasicGenome();
    
    // Store initial connections
    auto initialConnections = genome.get_phenotype()._orderedConnections;
    
    // Add connection between existing nodes
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    const auto& phenotype = genome.get_phenotype();
    
    // All original connections should still exist
    EXPECT_EQ(phenotype._orderedConnections.size(), initialConnections.size() + 1);
    
    // Verify original connections still exist (though indices might have changed)
    validateIndexMapping(genome);
}

// Edge Cases
TEST_F(PhenotypeUpdateConnectionTest, ConnectionToFromPreviouslyExcludedHidden_IncludesNode) {
    auto genome = createGenomeWithUnconnectedHidden();
    
    // Initially H4 should not be in phenotype
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    // Add connection involving H4
    addConnectionToGenome(genome, 13, 4, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    // Now H4 should be included
    EXPECT_EQ(genome.get_phenotype()._nodeGeneAttributes.size(), initialNodes + 1);
    validateNewConnectionExists(genome, 4, 3, 4.0f);
    validateIndexMapping(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, InputToOutputConnection_NoIntermediateNodes) {
    auto genome = createBasicGenome();
    
    // Add direct I1->O3 connection
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    size_t initialNodes = genome.get_phenotype()._nodeGeneAttributes.size();
    
    Operator::phenotypeUpdateConnection(genome);
    
    // Should not add new nodes
    EXPECT_EQ(genome.get_phenotype()._nodeGeneAttributes.size(), initialNodes);
    
    validateNewConnectionExists(genome, 1, 3, 4.0f);
    validatePhenotypeStructure(genome);
}

// State Consistency Tests
TEST_F(PhenotypeUpdateConnectionTest, DeltasCleared_AfterOperation) {
    auto genome = createBasicGenome();
    
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    validateDeltasEmpty(genome);
}

TEST_F(PhenotypeUpdateConnectionTest, ConnectionReferencesCorrectIndices) {
    auto genome = createBasicGenome();
    
    addConnectionToGenome(genome, 13, 1, 3, 4.0f, true);
    setConnectionDeltas(genome, {13});
    
    Operator::phenotypeUpdateConnection(genome);
    
    validatePhenotypeStructure(genome);
    validateIndexMapping(genome);
}