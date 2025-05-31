#include <gtest/gtest.h>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"

class GenomeTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}

    // Helper to create valid GenomeParams for basic testing
    GenomeParams createValidGenomeParams() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE},
            {ActivationType::SIGMOID}, 
            {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1};
        params._sourceNodeHistoryIDs = {1};
        params._targetNodeHistoryIDs = {2};
        params._connectionAttributes = {{1.5f, true}};
        return params;
    }

    // Helper to create valid RawGenomeParams for basic testing
    RawGenomeParams createValidRawGenomeParams() {
        static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
        static NodeGene hiddenNode(2, NodeType::HIDDEN, {ActivationType::SIGMOID});
        static NodeGene outputNode(3, NodeType::OUTPUT, {ActivationType::NONE});
        static ConnectionGene conn(1, inputNode, hiddenNode, {1.5f, true});
        
        RawGenomeParams params;
        params._nodeGenes = {&inputNode, &hiddenNode, &outputNode};
        params._rawConnectionGeneData = {conn.get_rawData()};
        return params;
    }
};

// ============================================================================
// GenomeParams Constructor Validation Tests
// ============================================================================

TEST_F(GenomeTest, GenomeParams_ValidConstruction) {
    GenomeParams params = createValidGenomeParams();
    EXPECT_NO_THROW(Genome genome(params));
}

TEST_F(GenomeTest, GenomeParams_NodeArraySizeMismatch_HistoryIDTypes) {
    GenomeParams params = createValidGenomeParams();
    params._nodeTypes.pop_back(); // Remove one element
    EXPECT_DEATH(Genome genome(params), "Node parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_NodeArraySizeMismatch_HistoryIDAttributes) {
    GenomeParams params = createValidGenomeParams();
    params._nodeAttributes.pop_back(); // Remove one element  
    EXPECT_DEATH(Genome genome(params), "Node parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_NodeArraySizeMismatch_TypesAttributes) {
    GenomeParams params = createValidGenomeParams();
    params._nodeHistoryIDs.pop_back(); // Remove one element
    EXPECT_DEATH(Genome genome(params), "Node parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_ConnectionArraySizeMismatch_HistoryIDSource) {
    GenomeParams params = createValidGenomeParams();
    params._sourceNodeHistoryIDs.pop_back(); // Remove one element
    EXPECT_DEATH(Genome genome(params), "Connection parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_ConnectionArraySizeMismatch_HistoryIDTarget) {
    GenomeParams params = createValidGenomeParams();
    params._targetNodeHistoryIDs.pop_back(); // Remove one element
    EXPECT_DEATH(Genome genome(params), "Connection parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_ConnectionArraySizeMismatch_HistoryIDAttributes) {
    GenomeParams params = createValidGenomeParams();
    params._connectionAttributes.pop_back(); // Remove one element
    EXPECT_DEATH(Genome genome(params), "Connection parameter arrays must be of equal size");
}

TEST_F(GenomeTest, GenomeParams_NonexistentSourceNode) {
    GenomeParams params = createValidGenomeParams();
    params._sourceNodeHistoryIDs[0] = 999; // Reference non-existent node
    EXPECT_DEATH(Genome genome(params), "Connection references nonexistent source node");
}

TEST_F(GenomeTest, GenomeParams_NonexistentTargetNode) {
    GenomeParams params = createValidGenomeParams();
    params._targetNodeHistoryIDs[0] = 999; // Reference non-existent node
    EXPECT_DEATH(Genome genome(params), "Connection references nonexistent target node");
}

TEST_F(GenomeTest, GenomeParams_OutputNodeAsSource) {
    GenomeParams params = createValidGenomeParams();
    // Try to use OUTPUT node (ID 3) as source
    params._sourceNodeHistoryIDs[0] = 3;
    params._targetNodeHistoryIDs[0] = 2;
    EXPECT_DEATH(Genome genome(params), "Output nodes cannot be connection sources");
}

TEST_F(GenomeTest, GenomeParams_InputNodeAsTarget) {
    GenomeParams params = createValidGenomeParams();
    // Try to use INPUT node (ID 1) as target
    params._sourceNodeHistoryIDs[0] = 2;
    params._targetNodeHistoryIDs[0] = 1;
    EXPECT_DEATH(Genome genome(params), "Input nodes cannot be connection targets");
}

TEST_F(GenomeTest, GenomeParams_DuplicateNodeHistoryIDs) {
    GenomeParams params = createValidGenomeParams();
    params._nodeHistoryIDs[1] = params._nodeHistoryIDs[0]; // Duplicate ID - node at index 1 now has same ID as node at index 0
    // Update connection to reference the duplicated ID (which still exists)
    params._targetNodeHistoryIDs[0] = params._nodeHistoryIDs[0]; // Use the duplicated ID
    // Note: Current implementation doesn't explicitly check for duplicates in GenomeParams constructor
    // This would create undefined behavior in the connection mapping, but last occurrence should win
    EXPECT_NO_THROW(Genome genome(params)); 
}

TEST_F(GenomeTest, GenomeParams_EmptyNodeVectors) {
    GenomeParams params;
    // All vectors empty - should succeed (minimal genome)
    EXPECT_NO_THROW(Genome genome(params));
}

TEST_F(GenomeTest, GenomeParams_NodesWithoutConnections) {
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2};
    params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    // No connections
    EXPECT_NO_THROW(Genome genome(params));
}

TEST_F(GenomeTest, GenomeParams_MultipleConnectionsSameNodes) {
    GenomeParams params = createValidGenomeParams();
    // Add another connection between same nodes
    params._connectionHistoryIDs.push_back(2);
    params._sourceNodeHistoryIDs.push_back(1);
    params._targetNodeHistoryIDs.push_back(2);
    params._connectionAttributes.push_back({2.0f, false});
    EXPECT_NO_THROW(Genome genome(params));
}

// ============================================================================
// RawGenomeParams Constructor Validation Tests  
// ============================================================================

TEST_F(GenomeTest, RawGenomeParams_ValidConstruction) {
    RawGenomeParams params = createValidRawGenomeParams();
    EXPECT_NO_THROW(Genome genome(params));
}

TEST_F(GenomeTest, RawGenomeParams_EmptyNodeGenes) {
    RawGenomeParams params;
    // Empty node genes vector
    EXPECT_DEATH(Genome genome(params), "Raw node gene data cannot be empty");
}

TEST_F(GenomeTest, RawGenomeParams_NullNodeGenePointer) {
    RawGenomeParams params = createValidRawGenomeParams();
    params._nodeGenes[0] = nullptr; // Null pointer
    EXPECT_DEATH(Genome genome(params), "Raw node data pointer cannot be null");
}

TEST_F(GenomeTest, RawGenomeParams_NullConnectionDataPointer) {
    RawGenomeParams params = createValidRawGenomeParams();
    params._rawConnectionGeneData[0] = nullptr; // Null pointer
    EXPECT_DEATH(Genome genome(params), "Raw connection data pointer cannot be null");
}

TEST_F(GenomeTest, RawGenomeParams_DuplicateNodeHistoryIDs) {
    static NodeGene node1(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene node2(1, NodeType::OUTPUT, {ActivationType::NONE}); // Same ID as node1
    
    RawGenomeParams params;
    params._nodeGenes = {&node1, &node2};
    EXPECT_DEATH(Genome genome(params), "Duplicate node history ID found in raw data");
}

TEST_F(GenomeTest, RawGenomeParams_ConnectionReferencesNonexistentNode) {
    static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene outputNode(3, NodeType::OUTPUT, {ActivationType::NONE});
    static NodeGene hiddenNode(2, NodeType::HIDDEN, {ActivationType::SIGMOID});
    // Create connection that references node ID 4 (doesn't exist)
    static ConnectionGene conn(1, hiddenNode, outputNode, {1.0f, true});
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &outputNode}; // Don't include hiddenNode (ID 2)
    params._rawConnectionGeneData = {conn.get_rawData()}; // Connection references ID 2
    EXPECT_DEATH(Genome genome(params), "Connection references nonexistent source node");
}

TEST_F(GenomeTest, RawGenomeParams_OutputNodeAsSource) {
    static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene outputNode(2, NodeType::OUTPUT, {ActivationType::NONE});
    static ConnectionGene conn(1, outputNode, inputNode, {1.0f, true}); // OUTPUT as source
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &outputNode};
    params._rawConnectionGeneData = {conn.get_rawData()};
    EXPECT_DEATH(Genome genome(params), "Output nodes cannot be connection sources");
}

TEST_F(GenomeTest, RawGenomeParams_InputNodeAsTarget) {
    static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene hiddenNode(2, NodeType::HIDDEN, {ActivationType::SIGMOID});
    static ConnectionGene conn(1, hiddenNode, inputNode, {1.0f, true}); // INPUT as target
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &hiddenNode};
    params._rawConnectionGeneData = {conn.get_rawData()};
    EXPECT_DEATH(Genome genome(params), "Input nodes cannot be connection targets");
}

// ============================================================================
// Constructor Equivalence Tests
// ============================================================================

TEST_F(GenomeTest, ConstructorEquivalence_BasicGenome) {
    // Create identical genome using both constructors
    GenomeParams genomeParams = createValidGenomeParams();
    RawGenomeParams rawParams = createValidRawGenomeParams();
    
    Genome genome1(genomeParams);
    Genome genome2(rawParams);
    
    // Verify same structure
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome2.get_connectionGenes().size());
    
    // Verify same node data
    for (size_t i = 0; i < genome1.get_nodeGenes().size(); ++i) {
        const auto& node1 = genome1.get_nodeGenes()[i];
        const auto& node2 = genome2.get_nodeGenes()[i];
        EXPECT_EQ(node1.get_historyID(), node2.get_historyID());
        EXPECT_EQ(node1.get_type(), node2.get_type());
        EXPECT_EQ(node1.get_attributes().activationType, node2.get_attributes().activationType);
    }
    
    // Verify same connection data
    for (size_t i = 0; i < genome1.get_connectionGenes().size(); ++i) {
        const auto& conn1 = genome1.get_connectionGenes()[i];
        const auto& conn2 = genome2.get_connectionGenes()[i];
        EXPECT_EQ(conn1.get_historyID(), conn2.get_historyID());
        EXPECT_EQ(conn1.get_sourceNodeGene().get_historyID(), conn2.get_sourceNodeGene().get_historyID());
        EXPECT_EQ(conn1.get_targetNodeGene().get_historyID(), conn2.get_targetNodeGene().get_historyID());
        EXPECT_EQ(conn1.get_attributes().weight, conn2.get_attributes().weight);
        EXPECT_EQ(conn1.get_attributes().enabled, conn2.get_attributes().enabled);
    }
}

TEST_F(GenomeTest, ConstructorEquivalence_ComplexGenome) {
    // Create a more complex genome for thorough testing
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2, 3, 4, 5};
    params._nodeTypes = {NodeType::INPUT, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
    params._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}, 
        {ActivationType::TANH}, {ActivationType::NONE}
    };
    params._connectionHistoryIDs = {1, 2, 3, 4};
    params._sourceNodeHistoryIDs = {1, 2, 3, 4};
    params._targetNodeHistoryIDs = {3, 4, 5, 5};
    params._connectionAttributes = {
        {1.0f, true}, {-0.5f, false}, {2.3f, true}, {0.8f, true}
    };
    
    // Create equivalent raw params
    static std::vector<NodeGene> staticNodes;
    staticNodes.clear();
    for (size_t i = 0; i < params._nodeHistoryIDs.size(); ++i) {
        staticNodes.emplace_back(params._nodeHistoryIDs[i], params._nodeTypes[i], params._nodeAttributes[i]);
    }
    
    static std::vector<ConnectionGene> staticConns;
    staticConns.clear();
    std::unordered_map<uint32_t, const NodeGene*> nodeMap;
    for (const auto& node : staticNodes) {
        nodeMap[node.get_historyID()] = &node;
    }
    
    for (size_t i = 0; i < params._connectionHistoryIDs.size(); ++i) {
        const NodeGene* source = nodeMap[params._sourceNodeHistoryIDs[i]];
        const NodeGene* target = nodeMap[params._targetNodeHistoryIDs[i]];
        staticConns.emplace_back(params._connectionHistoryIDs[i], *source, *target, params._connectionAttributes[i]);
    }
    
    RawGenomeParams rawParams;
    for (const auto& node : staticNodes) {
        rawParams._nodeGenes.push_back(&node);
    }
    for (const auto& conn : staticConns) {
        rawParams._rawConnectionGeneData.push_back(conn.get_rawData());
    }
    
    Genome genome1(params);
    Genome genome2(rawParams);
    
    // Verify equivalence
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome2.get_connectionGenes().size());
    
    for (size_t i = 0; i < genome1.get_nodeGenes().size(); ++i) {
        const auto& node1 = genome1.get_nodeGenes()[i];
        const auto& node2 = genome2.get_nodeGenes()[i];
        EXPECT_EQ(node1.get_historyID(), node2.get_historyID());
        EXPECT_EQ(node1.get_type(), node2.get_type());
        EXPECT_EQ(node1.get_attributes().activationType, node2.get_attributes().activationType);
    }
    
    for (size_t i = 0; i < genome1.get_connectionGenes().size(); ++i) {
        const auto& conn1 = genome1.get_connectionGenes()[i];
        const auto& conn2 = genome2.get_connectionGenes()[i];
        EXPECT_EQ(conn1.get_historyID(), conn2.get_historyID());
        EXPECT_EQ(conn1.get_sourceNodeGene().get_historyID(), conn2.get_sourceNodeGene().get_historyID());
        EXPECT_EQ(conn1.get_targetNodeGene().get_historyID(), conn2.get_targetNodeGene().get_historyID());
        EXPECT_FLOAT_EQ(conn1.get_attributes().weight, conn2.get_attributes().weight);
        EXPECT_EQ(conn1.get_attributes().enabled, conn2.get_attributes().enabled);
    }
}

// ============================================================================
// Copy Constructor Tests
// ============================================================================

TEST_F(GenomeTest, CopyConstructor_BasicCopy) {
    GenomeParams params = createValidGenomeParams();
    Genome original(params);
    Genome copy(original);
    
    // Verify same structure
    EXPECT_EQ(original.get_nodeGenes().size(), copy.get_nodeGenes().size());
    EXPECT_EQ(original.get_connectionGenes().size(), copy.get_connectionGenes().size());
    
    // Verify same data but different objects
    for (size_t i = 0; i < original.get_nodeGenes().size(); ++i) {
        const auto& origNode = original.get_nodeGenes()[i];
        const auto& copyNode = copy.get_nodeGenes()[i];
        
        // Same data
        EXPECT_EQ(origNode.get_historyID(), copyNode.get_historyID());
        EXPECT_EQ(origNode.get_type(), copyNode.get_type());
        EXPECT_EQ(origNode.get_attributes().activationType, copyNode.get_attributes().activationType);
        
        // Different objects
        EXPECT_NE(&origNode, &copyNode);
    }
    
    for (size_t i = 0; i < original.get_connectionGenes().size(); ++i) {
        const auto& origConn = original.get_connectionGenes()[i];
        const auto& copyConn = copy.get_connectionGenes()[i];
        
        // Same data
        EXPECT_EQ(origConn.get_historyID(), copyConn.get_historyID());
        EXPECT_EQ(origConn.get_attributes().weight, copyConn.get_attributes().weight);
        EXPECT_EQ(origConn.get_attributes().enabled, copyConn.get_attributes().enabled);
        
        // Different objects
        EXPECT_NE(&origConn, &copyConn);
        EXPECT_NE(&origConn.get_sourceNodeGene(), &copyConn.get_sourceNodeGene());
        EXPECT_NE(&origConn.get_targetNodeGene(), &copyConn.get_targetNodeGene());
        
        // But same history IDs
        EXPECT_EQ(origConn.get_sourceNodeGene().get_historyID(), copyConn.get_sourceNodeGene().get_historyID());
        EXPECT_EQ(origConn.get_targetNodeGene().get_historyID(), copyConn.get_targetNodeGene().get_historyID());
    }
}

TEST_F(GenomeTest, CopyConstructor_ReferenceRemapping) {
    // Create genome with multiple connections referencing same nodes
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2, 3};
    params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
    params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}};
    params._connectionHistoryIDs = {1, 2, 3};
    params._sourceNodeHistoryIDs = {1, 1, 2}; // Node 1 is source for multiple connections
    params._targetNodeHistoryIDs = {2, 3, 3}; // Node 3 is target for multiple connections
    params._connectionAttributes = {{1.0f, true}, {2.0f, true}, {3.0f, true}};
    
    Genome original(params);
    Genome copy(original);
    
    // Verify connections reference the correct copied nodes
    const auto& origConns = original.get_connectionGenes();
    const auto& copyConns = copy.get_connectionGenes();
    const auto& copyNodes = copy.get_nodeGenes();
    
    for (size_t i = 0; i < origConns.size(); ++i) {
        const auto& origConn = origConns[i];
        const auto& copyConn = copyConns[i];
        
        // Find the corresponding nodes in the copy by history ID
        const NodeGene* expectedSource = nullptr;
        const NodeGene* expectedTarget = nullptr;
        
        for (const auto& node : copyNodes) {
            if (node.get_historyID() == origConn.get_sourceNodeGene().get_historyID()) {
                expectedSource = &node;
            }
            if (node.get_historyID() == origConn.get_targetNodeGene().get_historyID()) {
                expectedTarget = &node;
            }
        }
        
        ASSERT_NE(expectedSource, nullptr);
        ASSERT_NE(expectedTarget, nullptr);
        
        // Verify the copied connection references the correct copied nodes
        EXPECT_EQ(&copyConn.get_sourceNodeGene(), expectedSource);
        EXPECT_EQ(&copyConn.get_targetNodeGene(), expectedTarget);
    }
}

TEST_F(GenomeTest, CopyConstructor_ModificationIndependence) {
    GenomeParams params = createValidGenomeParams();
    Genome original(params);
    Genome copy(original);
    
    // Modify original's mutable parts (would need access to mutable references)
    // Since we can't easily access mutable references in current API,
    // we verify that the objects are independent by checking addresses
    
    const auto& origNodes = original.get_nodeGenes();
    const auto& copyNodes = copy.get_nodeGenes();
    
    for (size_t i = 0; i < origNodes.size(); ++i) {
        EXPECT_NE(&origNodes[i], &copyNodes[i]) << "Node " << i << " should be independent";
    }
    
    const auto& origConns = original.get_connectionGenes();
    const auto& copyConns = copy.get_connectionGenes();
    
    for (size_t i = 0; i < origConns.size(); ++i) {
        EXPECT_NE(&origConns[i], &copyConns[i]) << "Connection " << i << " should be independent";
    }
}

TEST_F(GenomeTest, CopyConstructor_EmptyGenome) {
    GenomeParams params; // Empty genome
    Genome original(params);
    Genome copy(original);
    
    EXPECT_EQ(original.get_nodeGenes().size(), copy.get_nodeGenes().size());
    EXPECT_EQ(original.get_connectionGenes().size(), copy.get_connectionGenes().size());
    EXPECT_TRUE(copy.get_nodeGenes().empty());
    EXPECT_TRUE(copy.get_connectionGenes().empty());
}

TEST_F(GenomeTest, CopyConstructor_DisconnectedComponents) {
    // Create genome with nodes that aren't connected
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2, 3, 4};
    params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT, NodeType::BIAS};
    params._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::SIGMOID}, 
        {ActivationType::NONE}, {ActivationType::NONE}
    };
    params._connectionHistoryIDs = {1}; // Only one connection
    params._sourceNodeHistoryIDs = {1};
    params._targetNodeHistoryIDs = {3}; // Skip hidden node 2, ignore bias node 4
    params._connectionAttributes = {{1.0f, true}};
    
    Genome original(params);
    Genome copy(original);
    
    EXPECT_EQ(original.get_nodeGenes().size(), copy.get_nodeGenes().size());
    EXPECT_EQ(original.get_connectionGenes().size(), copy.get_connectionGenes().size());
    
    // Verify all nodes copied correctly, even disconnected ones
    for (size_t i = 0; i < original.get_nodeGenes().size(); ++i) {
        const auto& origNode = original.get_nodeGenes()[i];
        const auto& copyNode = copy.get_nodeGenes()[i];
        EXPECT_EQ(origNode.get_historyID(), copyNode.get_historyID());
        EXPECT_EQ(origNode.get_type(), copyNode.get_type());
        EXPECT_NE(&origNode, &copyNode);
    }
}

// ============================================================================
// Copy Assignment Tests
// ============================================================================

TEST_F(GenomeTest, CopyAssignment_BasicAssignment) {
    GenomeParams params1 = createValidGenomeParams();
    GenomeParams params2;
    params2._nodeHistoryIDs = {10, 20};
    params2._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params2._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    // No connections in params2
    
    Genome genome1(params1);
    Genome genome2(params2);
    
    // Verify initial difference
    EXPECT_NE(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    
    // Perform assignment
    genome2 = genome1;
    
    // Verify assignment worked
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome2.get_connectionGenes().size());
    
    // Verify data copied correctly
    for (size_t i = 0; i < genome1.get_nodeGenes().size(); ++i) {
        const auto& node1 = genome1.get_nodeGenes()[i];
        const auto& node2 = genome2.get_nodeGenes()[i];
        EXPECT_EQ(node1.get_historyID(), node2.get_historyID());
        EXPECT_EQ(node1.get_type(), node2.get_type());
        EXPECT_EQ(node1.get_attributes().activationType, node2.get_attributes().activationType);
        EXPECT_NE(&node1, &node2); // Different objects
    }
}

TEST_F(GenomeTest, CopyAssignment_SelfAssignment) {
    GenomeParams params = createValidGenomeParams();
    Genome genome(params);
    
    // Store original addresses to verify no change
    const void* origNodePtr = &genome.get_nodeGenes()[0];
    const void* origConnPtr = &genome.get_connectionGenes()[0];
    size_t origNodeCount = genome.get_nodeGenes().size();
    size_t origConnCount = genome.get_connectionGenes().size();
    
    // Self-assignment
    genome = genome;
    
    // Verify no change in structure
    EXPECT_EQ(origNodeCount, genome.get_nodeGenes().size());
    EXPECT_EQ(origConnCount, genome.get_connectionGenes().size());
    
    // Verify same objects (no unnecessary copying)
    EXPECT_EQ(origNodePtr, &genome.get_nodeGenes()[0]);
    EXPECT_EQ(origConnPtr, &genome.get_connectionGenes()[0]);
}

TEST_F(GenomeTest, CopyAssignment_PhenotypeClearing) {
    GenomeParams params = createValidGenomeParams();
    Genome genome1(params);
    Genome genome2(params);
    
    // Force phenotype construction on genome2
    genome2.constructPhenotype();
    EXPECT_NE(genome2.get_phenotype(), nullptr);
    
    // Assignment should clear phenotype
    genome2 = genome1;
    EXPECT_EQ(genome2.get_phenotype(), nullptr);
}

TEST_F(GenomeTest, CopyAssignment_StateCleanup) {
    // Create large genome
    GenomeParams largeParams;
    for (uint32_t i = 1; i <= 100; ++i) {
        largeParams._nodeHistoryIDs.push_back(i);
        largeParams._nodeTypes.push_back(i <= 10 ? NodeType::INPUT : 
                                       i > 90 ? NodeType::OUTPUT : NodeType::HIDDEN);
        largeParams._nodeAttributes.push_back({ActivationType::SIGMOID});
    }
    
    GenomeParams smallParams = createValidGenomeParams();
    
    Genome largeGenome(largeParams);
    Genome smallGenome(smallParams);
    
    // Assign small to large (should properly clean up large state)
    largeGenome = smallGenome;
    
    EXPECT_EQ(largeGenome.get_nodeGenes().size(), smallGenome.get_nodeGenes().size());
    EXPECT_EQ(largeGenome.get_connectionGenes().size(), smallGenome.get_connectionGenes().size());
}

TEST_F(GenomeTest, CopyAssignment_AssignmentChaining) {
    GenomeParams params1 = createValidGenomeParams();
    
    GenomeParams params2;
    params2._nodeHistoryIDs = {10, 20};
    params2._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params2._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    
    GenomeParams params3;
    params3._nodeHistoryIDs = {100};
    params3._nodeTypes = {NodeType::BIAS};
    params3._nodeAttributes = {{ActivationType::NONE}};
    
    Genome genome1(params1);
    Genome genome2(params2);
    Genome genome3(params3);
    
    // Chain assignment
    genome3 = genome2 = genome1;
    
    // All should be equal to genome1
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome3.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome2.get_connectionGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome3.get_connectionGenes().size());
    
    // Verify data consistency
    for (size_t i = 0; i < genome1.get_nodeGenes().size(); ++i) {
        EXPECT_EQ(genome1.get_nodeGenes()[i].get_historyID(), genome2.get_nodeGenes()[i].get_historyID());
        EXPECT_EQ(genome1.get_nodeGenes()[i].get_historyID(), genome3.get_nodeGenes()[i].get_historyID());
    }
}

TEST_F(GenomeTest, CopyAssignment_ComplexReferenceRemapping) {
    // Create complex genome with multiple connection patterns
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2, 3, 4, 5};
    params._nodeTypes = {NodeType::INPUT, NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
    params._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID},
        {ActivationType::TANH}, {ActivationType::NONE}
    };
    params._connectionHistoryIDs = {1, 2, 3, 4, 5, 6};
    params._sourceNodeHistoryIDs = {1, 2, 1, 3, 4, 3};
    params._targetNodeHistoryIDs = {3, 3, 4, 4, 5, 5};
    params._connectionAttributes = {
        {1.0f, true}, {2.0f, true}, {3.0f, false}, {4.0f, true}, {5.0f, true}, {6.0f, false}
    };
    
    GenomeParams simpleParams = createValidGenomeParams();
    
    Genome complexGenome(params);
    Genome simpleGenome(simpleParams);
    
    // Assign complex to simple
    simpleGenome = complexGenome;
    
    // Verify reference remapping worked correctly
    const auto& nodes = simpleGenome.get_nodeGenes();
    const auto& conns = simpleGenome.get_connectionGenes();
    
    for (const auto& conn : conns) {
        // Find source and target nodes by history ID
        const NodeGene* sourceFound = nullptr;
        const NodeGene* targetFound = nullptr;
        
        for (const auto& node : nodes) {
            if (node.get_historyID() == conn.get_sourceNodeGene().get_historyID()) {
                sourceFound = &node;
            }
            if (node.get_historyID() == conn.get_targetNodeGene().get_historyID()) {
                targetFound = &node;
            }
        }
        
        ASSERT_NE(sourceFound, nullptr) << "Source node not found for connection " << conn.get_historyID();
        ASSERT_NE(targetFound, nullptr) << "Target node not found for connection " << conn.get_historyID();
        
        // Verify connection references the correct nodes from our copied set
        EXPECT_EQ(&conn.get_sourceNodeGene(), sourceFound);
        EXPECT_EQ(&conn.get_targetNodeGene(), targetFound);
    }
}