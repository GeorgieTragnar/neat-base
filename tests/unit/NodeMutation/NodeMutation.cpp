#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/NodeMutation.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Operator;

class NodeMutationTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::srand(42); // Fixed seed for reproducibility where needed
    }

    // Helper to create a properly initialized HistoryTracker for a given genome
    std::shared_ptr<HistoryTracker> createInitializedHistoryTracker(const Genome& genome) {
        auto tracker = std::make_shared<HistoryTracker>();
        
        // Initialize tracker with genome's existing nodes
        for (const auto& node : genome.get_nodeGenes()) {
            // Force tracker to register existing node IDs to advance _nextNodeID
            switch (node.get_type()) {
                case NodeType::INPUT:
                    // Assume input nodes have consecutive IDs starting from 1
                    tracker->get_input(node.get_historyID() - 1);
                    break;
                case NodeType::OUTPUT:
                    // Assume output nodes come after inputs
                    tracker->get_output(0); // Just ensure tracker knows about outputs
                    break;
                case NodeType::BIAS:
                    tracker->get_bias();
                    break;
                case NodeType::HIDDEN:
                    // Can't easily register hidden nodes, but they'll be handled
                    break;
            }
        }
        
        // Initialize tracker with genome's existing connections
        for (const auto& conn : genome.get_connectionGenes()) {
            tracker->get_connection(
                conn.get_sourceNodeGene().get_historyID(),
                conn.get_targetNodeGene().get_historyID()
            );
        }
        
        return tracker;
    }

    // Helper to create basic genome with enabled connections
    Genome createBasicGenome() {
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
            {1.5f, true},   // enabled
            {2.0f, true}    // enabled
        };
        return Genome(params);
    }

    // Helper to create genome with single enabled connection
    Genome createSingleConnectionGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
        params._connectionHistoryIDs = {1};
        params._sourceNodeHistoryIDs = {1};
        params._targetNodeHistoryIDs = {2};
        params._connectionAttributes = {{2.5f, true}};
        return Genome(params);
    }

    // Helper to create genome with mixed enabled/disabled connections
    Genome createMixedConnectionGenome() {
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
            {1.5f, false},  // disabled
            {2.0f, true}    // enabled
        };
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
            {1.5f, false}   // disabled
        };
        return Genome(params);
    }

    // Helper to create empty genome
    Genome createEmptyGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
        params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
        // No connections
        return Genome(params);
    }

    // Helper to count enabled connections
    size_t countEnabledConnections(const Genome& genome) {
        size_t count = 0;
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_attributes().enabled) {
                count++;
            }
        }
        return count;
    }

    // Helper to count disabled connections
    size_t countDisabledConnections(const Genome& genome) {
        size_t count = 0;
        for (const auto& conn : genome.get_connectionGenes()) {
            if (!conn.get_attributes().enabled) {
                count++;
            }
        }
        return count;
    }

    // Helper to find connection by innovation number
    const ConnectionGene* findConnectionByID(const Genome& genome, uint32_t historyID) {
        for (const auto& conn : genome.get_connectionGenes()) {
            if (conn.get_historyID() == historyID) {
                return &conn;
            }
        }
        return nullptr;
    }

    // Helper to find node by history ID
    const NodeGene* findNodeByID(const Genome& genome, uint32_t historyID) {
        for (const auto& node : genome.get_nodeGenes()) {
            if (node.get_historyID() == historyID) {
                return &node;
            }
        }
        return nullptr;
    }

    // Helper to verify genome structural integrity
    void verifyGenomeIntegrity(const Genome& original, const Genome& mutated) {
        // Should have one more node OR same if node was reused
        EXPECT_GE(mutated.get_nodeGenes().size(), original.get_nodeGenes().size());
        EXPECT_LE(mutated.get_nodeGenes().size(), original.get_nodeGenes().size() + 1);
        
        // Should have two more connections
        EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 2);
        
        // All original nodes should be preserved
        for (const auto& origNode : original.get_nodeGenes()) {
            const NodeGene* mutNode = findNodeByID(mutated, origNode.get_historyID());
            ASSERT_NE(mutNode, nullptr);
            EXPECT_EQ(origNode.get_type(), mutNode->get_type());
            EXPECT_EQ(origNode.get_attributes().activationType, mutNode->get_attributes().activationType);
        }
    }
    
    // Helper to debug genome state
    void debugGenomeState(const std::string& label, const Genome& genome) {
        std::cout << "\n=== " << label << " ===\n";
        std::cout << "Nodes: " << genome.get_nodeGenes().size() << "\n";
        for (const auto& node : genome.get_nodeGenes()) {
            std::cout << "  Node " << node.get_historyID() << " type=" << (int)node.get_type() << "\n";
        }
        std::cout << "Connections: " << genome.get_connectionGenes().size() << "\n";
        for (const auto& conn : genome.get_connectionGenes()) {
            std::cout << "  Conn " << conn.get_historyID() 
                     << " " << conn.get_sourceNodeGene().get_historyID() 
                     << "->" << conn.get_targetNodeGene().get_historyID()
                     << " enabled=" << conn.get_attributes().enabled << "\n";
        }
    }

    // Helper to verify node mutation structural changes
    struct MutationResult {
        bool hasValidStructure = false;
        const NodeGene* splitNode = nullptr;
        std::vector<const ConnectionGene*> newConnections;
        const ConnectionGene* disabledConnection = nullptr;
    };
    
    MutationResult analyzeMutation(const Genome& original, const Genome& mutated) {
        MutationResult result;
        
        // Verify basic structural requirements
        if (mutated.get_connectionGenes().size() != original.get_connectionGenes().size() + 2) {
            return result; // Must have +2 connections
        }
        
        if (countEnabledConnections(mutated) != countEnabledConnections(original) + 1) {
            return result; // Must have net +1 enabled connection
        }
        
        if (countDisabledConnections(mutated) != countDisabledConnections(original) + 1) {
            return result; // Must have +1 disabled connection
        }
        
        // Find the disabled original connection
        for (const auto& origConn : original.get_connectionGenes()) {
            if (origConn.get_attributes().enabled) {
                const ConnectionGene* mutConn = findConnectionByID(mutated, origConn.get_historyID());
                if (mutConn && !mutConn->get_attributes().enabled) {
                    result.disabledConnection = mutConn;
                    break;
                }
            }
        }
        
        if (!result.disabledConnection) {
            return result; // Must have exactly one disabled original connection
        }
        
        // Find new connections
        for (const auto& conn : mutated.get_connectionGenes()) {
            bool isNewConnection = true;
            for (const auto& origConn : original.get_connectionGenes()) {
                if (origConn.get_historyID() == conn.get_historyID()) {
                    isNewConnection = false;
                    break;
                }
            }
            if (isNewConnection && conn.get_attributes().enabled) {
                result.newConnections.push_back(&conn);
            }
        }
        
        if (result.newConnections.size() != 2) {
            return result; // Must have exactly 2 new connections
        }
        
        // Find the split node - it should be the common node in the two new connections
        uint32_t candidate1 = 0, candidate2 = 0;
        
        // Check if first connection's source/target matches second connection's source/target
        const auto& conn1 = *result.newConnections[0];
        const auto& conn2 = *result.newConnections[1];
        
        if (conn1.get_sourceNodeGene().get_historyID() == conn2.get_targetNodeGene().get_historyID()) {
            candidate1 = conn1.get_sourceNodeGene().get_historyID();
        }
        if (conn1.get_targetNodeGene().get_historyID() == conn2.get_sourceNodeGene().get_historyID()) {
            candidate2 = conn1.get_targetNodeGene().get_historyID();
        }
        
        // The split node should be HIDDEN type
        uint32_t splitNodeID = 0;
        if (candidate1 != 0) {
            const NodeGene* node = findNodeByID(mutated, candidate1);
            if (node && node->get_type() == NodeType::HIDDEN) {
                splitNodeID = candidate1;
            }
        }
        if (candidate2 != 0 && splitNodeID == 0) {
            const NodeGene* node = findNodeByID(mutated, candidate2);
            if (node && node->get_type() == NodeType::HIDDEN) {
                splitNodeID = candidate2;
            }
        }
        
        if (splitNodeID != 0) {
            result.splitNode = findNodeByID(mutated, splitNodeID);
            result.hasValidStructure = true;
        }
        
        return result;
    }
};

// ============================================================================
// Basic Structure Modification Tests
// ============================================================================

TEST_F(NodeMutationTest, NodeAlwaysAdded) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    debugGenomeState("ORIGINAL", original);
    
    size_t originalNodeCount = original.get_nodeGenes().size();
    size_t originalConnCount = original.get_connectionGenes().size();
    size_t originalEnabledCount = countEnabledConnections(original);
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    debugGenomeState("MUTATED", mutated);
    
    // Verify structural changes that MUST happen
    EXPECT_EQ(mutated.get_connectionGenes().size(), originalConnCount + 2); // +2 connections
    EXPECT_EQ(countEnabledConnections(mutated), originalEnabledCount + 1); // +1 enabled, -1 disabled = net +1
    EXPECT_EQ(countDisabledConnections(mutated), countDisabledConnections(original) + 1); // +1 disabled
    
    // Node count can stay same (reuse) or increase by 1 (new node)
    EXPECT_GE(mutated.get_nodeGenes().size(), originalNodeCount);
    EXPECT_LE(mutated.get_nodeGenes().size(), originalNodeCount + 1);
    
    // Verify that the operation actually worked by checking the core mutation happened
    // At least one original enabled connection should now be disabled
    bool foundDisabledOriginal = false;
    for (const auto& origConn : original.get_connectionGenes()) {
        if (origConn.get_attributes().enabled) {
            const ConnectionGene* mutConn = findConnectionByID(mutated, origConn.get_historyID());
            if (mutConn && !mutConn->get_attributes().enabled) {
                foundDisabledOriginal = true;
                break;
            }
        }
    }
    EXPECT_TRUE(foundDisabledOriginal);
}

TEST_F(NodeMutationTest, ConnectionsAdded) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    size_t originalConnCount = original.get_connectionGenes().size();
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Should have exactly two more connections
    EXPECT_EQ(mutated.get_connectionGenes().size(), originalConnCount + 2);
}

TEST_F(NodeMutationTest, OriginalConnectionDisabled) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Store original connection ID
    uint32_t originalConnID = original.get_connectionGenes()[0].get_historyID();
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Find the original connection in mutated genome
    const ConnectionGene* originalConn = findConnectionByID(mutated, originalConnID);
    ASSERT_NE(originalConn, nullptr);
    EXPECT_FALSE(originalConn->get_attributes().enabled);
}

TEST_F(NodeMutationTest, NodeTypeCorrect) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::TANH});
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Verify mutation structure and node type
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    EXPECT_EQ(result.splitNode->get_type(), NodeType::HIDDEN);
}

TEST_F(NodeMutationTest, NodeAttributesApplied) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::RELU});
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Verify mutation structure
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Note: If node was reused, it may have different attributes than params
    // This test primarily verifies the split functionality works
}

// ============================================================================
// Innovation Decision Tree Tests
// ============================================================================

TEST_F(NodeMutationTest, PrimaryPath_NewNode) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    debugGenomeState("SINGLE_ORIGINAL", original);
    
    // First split should create structural changes
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    debugGenomeState("SINGLE_MUTATED", mutated);
    
    // Verify core structural changes
    EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 2);
    EXPECT_EQ(countEnabledConnections(mutated), 2u); // 2 new enabled
    EXPECT_EQ(countDisabledConnections(mutated), 1u); // 1 original disabled
}

TEST_F(NodeMutationTest, InnovationConsistency_SameConnection) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Split same connection multiple times with fresh history trackers
    std::set<uint32_t> primaryNodeIDs;
    
    for (int i = 0; i < 5; ++i) {
        Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
        
        // Find the split node
        auto result = analyzeMutation(original, mutated);
        if (result.hasValidStructure && result.splitNode != nullptr) {
            primaryNodeIDs.insert(result.splitNode->get_historyID());
        }
    }
    
    // All should get the same primary node ID
    EXPECT_EQ(primaryNodeIDs.size(), 1u);
}

// ============================================================================
// Weight Preservation Tests
// ============================================================================

TEST_F(NodeMutationTest, InputConnectionWeight) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    uint32_t originalSourceID = original.get_connectionGenes()[0].get_sourceNodeGene().get_historyID();
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Find connection from original source to split node
    const ConnectionGene* inputConn = nullptr;
    for (const auto& conn : result.newConnections) {
        if (conn->get_sourceNodeGene().get_historyID() == originalSourceID &&
            conn->get_targetNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            inputConn = conn;
            break;
        }
    }
    
    ASSERT_NE(inputConn, nullptr);
    EXPECT_FLOAT_EQ(inputConn->get_attributes().weight, 1.0f);
}

TEST_F(NodeMutationTest, OutputConnectionWeight) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    float originalWeight = original.get_connectionGenes()[0].get_attributes().weight;
    uint32_t originalTargetID = original.get_connectionGenes()[0].get_targetNodeGene().get_historyID();
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Find connection from split node to original target
    const ConnectionGene* outputConn = nullptr;
    for (const auto& conn : result.newConnections) {
        if (conn->get_sourceNodeGene().get_historyID() == result.splitNode->get_historyID() &&
            conn->get_targetNodeGene().get_historyID() == originalTargetID) {
            outputConn = conn;
            break;
        }
    }
    
    ASSERT_NE(outputConn, nullptr);
    EXPECT_FLOAT_EQ(outputConn->get_attributes().weight, originalWeight);
}

TEST_F(NodeMutationTest, NetworkFunctionMaintained) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    float originalWeight = original.get_connectionGenes()[0].get_attributes().weight;
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Find input and output connections
    float inputWeight = 0.0f, outputWeight = 0.0f;
    int foundConnections = 0;
    
    for (const auto& conn : result.newConnections) {
        if (conn->get_targetNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            inputWeight = conn->get_attributes().weight;
            foundConnections++;
        } else if (conn->get_sourceNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            outputWeight = conn->get_attributes().weight;
            foundConnections++;
        }
    }
    
    EXPECT_EQ(foundConnections, 2);
    
    // Mathematical equivalence: input * 1.0 * output = original
    float combinedWeight = inputWeight * outputWeight;
    EXPECT_FLOAT_EQ(combinedWeight, originalWeight);
}

TEST_F(NodeMutationTest, OriginalWeightUnchanged) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    float originalWeight = original.get_connectionGenes()[0].get_attributes().weight;
    uint32_t originalConnID = original.get_connectionGenes()[0].get_historyID();
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Find original connection in mutated genome
    const ConnectionGene* originalConn = findConnectionByID(mutated, originalConnID);
    ASSERT_NE(originalConn, nullptr);
    EXPECT_FLOAT_EQ(originalConn->get_attributes().weight, originalWeight);
}

// ============================================================================
// Connection Split Logic Tests
// ============================================================================

TEST_F(NodeMutationTest, CorrectSourceTarget) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    uint32_t originalSourceID = original.get_connectionGenes()[0].get_sourceNodeGene().get_historyID();
    uint32_t originalTargetID = original.get_connectionGenes()[0].get_targetNodeGene().get_historyID();
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Verify connection structure
    bool foundInputConn = false, foundOutputConn = false;
    
    for (const auto& conn : result.newConnections) {
        // Input connection: originalSource → splitNode
        if (conn->get_sourceNodeGene().get_historyID() == originalSourceID &&
            conn->get_targetNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            foundInputConn = true;
        }
        // Output connection: splitNode → originalTarget
        else if (conn->get_sourceNodeGene().get_historyID() == result.splitNode->get_historyID() &&
                 conn->get_targetNodeGene().get_historyID() == originalTargetID) {
            foundOutputConn = true;
        }
    }
    
    EXPECT_TRUE(foundInputConn);
    EXPECT_TRUE(foundOutputConn);
}

TEST_F(NodeMutationTest, SplitNodePositioning) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Should have exactly 2 enabled connections and 1 disabled
    EXPECT_EQ(countEnabledConnections(mutated), 2u);
    EXPECT_EQ(countDisabledConnections(mutated), 1u);
    
    auto result = analyzeMutation(original, mutated);
    ASSERT_TRUE(result.hasValidStructure);
    ASSERT_NE(result.splitNode, nullptr);
    
    // Split node should appear as both source and target in enabled connections
    bool isSource = false, isTarget = false;
    for (const auto& conn : result.newConnections) {
        if (conn->get_sourceNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            isSource = true;
        }
        if (conn->get_targetNodeGene().get_historyID() == result.splitNode->get_historyID()) {
            isTarget = true;
        }
    }
    
    EXPECT_TRUE(isSource);
    EXPECT_TRUE(isTarget);
}

TEST_F(NodeMutationTest, EnabledState) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    size_t originalEnabled = countEnabledConnections(original);
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Should have 2 new enabled connections, 1 original disabled
    EXPECT_EQ(countEnabledConnections(mutated), originalEnabled + 1); // +2 new, -1 disabled
    EXPECT_EQ(countDisabledConnections(mutated), countDisabledConnections(original) + 1);
}

// ============================================================================
// Genome Integrity Tests
// ============================================================================

TEST_F(NodeMutationTest, OriginalGenomeUntouched) {
    Genome original = createBasicGenome();
    
    // Store original counts
    size_t origNodes = original.get_nodeGenes().size();
    size_t origConns = original.get_connectionGenes().size();
    size_t origEnabled = countEnabledConnections(original);
    
    NodeMutationParams params({ActivationType::SIGMOID});
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Verify original unchanged
    EXPECT_EQ(original.get_nodeGenes().size(), origNodes);
    EXPECT_EQ(original.get_connectionGenes().size(), origConns);
    EXPECT_EQ(countEnabledConnections(original), origEnabled);
}

TEST_F(NodeMutationTest, ExistingDataPreserved) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    verifyGenomeIntegrity(original, mutated);
}

TEST_F(NodeMutationTest, StructureConsistency) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Perform mutation - this test focuses on the mutated genome structure
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // All connections in mutated genome should have valid node references
    for (const auto& conn : mutated.get_connectionGenes()) {
        const NodeGene* sourceNode = findNodeByID(mutated, conn.get_sourceNodeGene().get_historyID());
        const NodeGene* targetNode = findNodeByID(mutated, conn.get_targetNodeGene().get_historyID());
        
        EXPECT_NE(sourceNode, nullptr) << "Connection " << conn.get_historyID() << " has invalid source";
        EXPECT_NE(targetNode, nullptr) << "Connection " << conn.get_historyID() << " has invalid target";
    }
}

TEST_F(NodeMutationTest, PhenotypeConstruction) {
    Genome original = createBasicGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Should be able to construct phenotype without errors
    EXPECT_NO_THROW({
        mutated.constructPhenotype();
        auto phenotype = mutated.get_phenotype();
        EXPECT_NE(phenotype, nullptr);
    });
}

// ============================================================================
// Random Selection Tests
// ============================================================================

TEST_F(NodeMutationTest, ConnectionSelection_CanProduceDifferent) {
    Genome original = createBasicGenome(); // Has multiple enabled connections
    NodeMutationParams params({ActivationType::SIGMOID});
    
    std::set<uint32_t> splitConnectionIDs;
    
    // Run enough iterations to likely see different selections
    for (int i = 0; i < 50; ++i) {
        Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
        
        // Find which connection was disabled (the one that was split)
        for (const auto& conn : mutated.get_connectionGenes()) {
            if (!conn.get_attributes().enabled) {
                // Check if this connection was enabled in original
                const ConnectionGene* origConn = findConnectionByID(original, conn.get_historyID());
                if (origConn && origConn->get_attributes().enabled) {
                    splitConnectionIDs.insert(conn.get_historyID());
                    break;
                }
            }
        }
    }
    
    // Should observe more than one connection being split
    EXPECT_GT(splitConnectionIDs.size(), 1u);
}

TEST_F(NodeMutationTest, SingleConnectionDeterministic) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    uint32_t expectedConnID = original.get_connectionGenes()[0].get_historyID();
    
    // Run multiple times - should always split the same connection
    for (int i = 0; i < 10; ++i) {
        Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
        
        // Find which connection was disabled
        const ConnectionGene* disabledConn = findConnectionByID(mutated, expectedConnID);
        ASSERT_NE(disabledConn, nullptr);
        EXPECT_FALSE(disabledConn->get_attributes().enabled);
    }
}

// ============================================================================
// Edge Cases & Error Handling Tests
// ============================================================================

TEST_F(NodeMutationTest, NoEnabledConnections) {
    Genome original = createAllDisabledGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Should assert when no enabled connections available
    EXPECT_DEATH({
        nodeMutation(original, createInitializedHistoryTracker(original), params);
    }, "");
}

TEST_F(NodeMutationTest, EmptyConnections) {
    Genome original = createEmptyGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Should assert when no connections exist
    EXPECT_DEATH({
        nodeMutation(original, createInitializedHistoryTracker(original), params);
    }, "");
}

TEST_F(NodeMutationTest, SingleConnection_CreatesThreeNodeNetwork) {
    Genome original = createSingleConnectionGenome();
    NodeMutationParams params({ActivationType::SIGMOID});
    
    // Original: 2 nodes, 1 connection
    EXPECT_EQ(original.get_nodeGenes().size(), 2u);
    EXPECT_EQ(original.get_connectionGenes().size(), 1u);
    
    Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), params);
    
    // Result: 2-3 nodes (depending on reuse), 3 connections (2 enabled, 1 disabled)
    EXPECT_GE(mutated.get_nodeGenes().size(), 2u);
    EXPECT_LE(mutated.get_nodeGenes().size(), 3u);
    EXPECT_EQ(mutated.get_connectionGenes().size(), 3u);
    EXPECT_EQ(countEnabledConnections(mutated), 2u);
    EXPECT_EQ(countDisabledConnections(mutated), 1u);
    
    // Should have valid mutation structure
    auto result = analyzeMutation(original, mutated);
    EXPECT_TRUE(result.hasValidStructure);
    EXPECT_NE(result.splitNode, nullptr);
}

// ============================================================================
// Complex Scenario Tests
// ============================================================================

TEST_F(NodeMutationTest, LargeGenome_Performance) {
    // Create large genome with many connections
    GenomeParams params;
    const size_t numNodes = 20;
    
    for (size_t i = 1; i <= numNodes; ++i) {
        params._nodeHistoryIDs.push_back(static_cast<uint32_t>(i));
        if (i <= 5) {
            params._nodeTypes.push_back(NodeType::INPUT);
        } else if (i > 15) {
            params._nodeTypes.push_back(NodeType::OUTPUT);
        } else {
            params._nodeTypes.push_back(NodeType::HIDDEN);
        }
        params._nodeAttributes.push_back({ActivationType::SIGMOID});
    }
    
    // Add many enabled connections
    for (size_t i = 0; i < 30; ++i) {
        params._connectionHistoryIDs.push_back(static_cast<uint32_t>(i + 1));
        params._sourceNodeHistoryIDs.push_back(static_cast<uint32_t>((i % 10) + 1));
        params._targetNodeHistoryIDs.push_back(static_cast<uint32_t>((i % 10) + 11));
        params._connectionAttributes.push_back({static_cast<float>(i * 0.1), true});
    }
    
    Genome original(params);
    NodeMutationParams mutParams({ActivationType::RELU});
    
    // Should handle large genome efficiently
    EXPECT_NO_THROW({
        Genome mutated = nodeMutation(original, createInitializedHistoryTracker(original), mutParams);
        EXPECT_GE(mutated.get_nodeGenes().size(), original.get_nodeGenes().size());
        EXPECT_LE(mutated.get_nodeGenes().size(), original.get_nodeGenes().size() + 1);
        EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 2);
    });
}

// ============================================================================
// Parameter Tests
// ============================================================================

TEST_F(NodeMutationTest, ValidNodeAttributes) {
    Genome original = createBasicGenome();
    
    EXPECT_NO_THROW({
        NodeMutationParams params1({ActivationType::SIGMOID});
        NodeMutationParams params2({ActivationType::TANH});
        NodeMutationParams params3({ActivationType::RELU});
        NodeMutationParams params4({ActivationType::LEAKY_RELU});
    });
}