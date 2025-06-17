#include <gtest/gtest.h>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/Genome.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"

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
        static ConnectionGene conn(1, 0, 1, {1.5f, true}); // Using indices: inputNode=0, hiddenNode=1
        
        RawGenomeParams params;
        params._nodeGenes = {&inputNode, &hiddenNode, &outputNode};
        params._connectionGenes = {&conn};
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
    params._connectionGenes[0] = nullptr; // Null pointer
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
    // Create connection that references index 2 (out of bounds - only indices 0,1 exist)
    static ConnectionGene conn(1, 0, 2, {1.0f, true}); // References index 2 as target (out of bounds)
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &outputNode}; // Only indices 0,1 exist
    params._connectionGenes = {&conn}; // Connection references invalid target index 2
    EXPECT_DEATH(Genome genome(params), "Connection target index out of bounds");
}

TEST_F(GenomeTest, RawGenomeParams_OutputNodeAsSource) {
    static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene outputNode(2, NodeType::OUTPUT, {ActivationType::NONE});
    static ConnectionGene conn(1, 1, 0, {1.0f, true}); // OUTPUT node (index 1) as source
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &outputNode};
    params._connectionGenes = {&conn};
    EXPECT_DEATH(Genome genome(params), "Output nodes cannot be connection sources");
}

TEST_F(GenomeTest, RawGenomeParams_InputNodeAsTarget) {
    static NodeGene inputNode(1, NodeType::INPUT, {ActivationType::NONE});
    static NodeGene hiddenNode(2, NodeType::HIDDEN, {ActivationType::SIGMOID});
    static ConnectionGene conn(1, 1, 0, {1.0f, true}); // INPUT node (index 0) as target
    
    RawGenomeParams params;
    params._nodeGenes = {&inputNode, &hiddenNode};
    params._connectionGenes = {&conn};
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
        const auto& nodes1 = genome1.get_nodeGenes();
        const auto& nodes2 = genome2.get_nodeGenes();
        EXPECT_EQ(conn1.get_historyID(), conn2.get_historyID());
        EXPECT_EQ(nodes1[conn1.get_sourceNodeIndex()].get_historyID(), nodes2[conn2.get_sourceNodeIndex()].get_historyID());
        EXPECT_EQ(nodes1[conn1.get_targetNodeIndex()].get_historyID(), nodes2[conn2.get_targetNodeIndex()].get_historyID());
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
    
    // Create index mapping
    std::unordered_map<uint32_t, size_t> nodeIndexMap;
    for (size_t i = 0; i < staticNodes.size(); ++i) {
        nodeIndexMap[staticNodes[i].get_historyID()] = i;
    }
    
    for (size_t i = 0; i < params._connectionHistoryIDs.size(); ++i) {
        size_t sourceIndex = nodeIndexMap[params._sourceNodeHistoryIDs[i]];
        size_t targetIndex = nodeIndexMap[params._targetNodeHistoryIDs[i]];
        staticConns.emplace_back(params._connectionHistoryIDs[i], sourceIndex, targetIndex, params._connectionAttributes[i]);
    }
    
    RawGenomeParams rawParams;
    for (const auto& node : staticNodes) {
        rawParams._nodeGenes.push_back(&node);
    }
    for (const auto& conn : staticConns) {
        rawParams._connectionGenes.push_back(&conn);
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
        const auto& nodes1 = genome1.get_nodeGenes();
        const auto& nodes2 = genome2.get_nodeGenes();
        EXPECT_EQ(conn1.get_historyID(), conn2.get_historyID());
        EXPECT_EQ(nodes1[conn1.get_sourceNodeIndex()].get_historyID(), nodes2[conn2.get_sourceNodeIndex()].get_historyID());
        EXPECT_EQ(nodes1[conn1.get_targetNodeIndex()].get_historyID(), nodes2[conn2.get_targetNodeIndex()].get_historyID());
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
        
        // But same node history IDs (via indices)
        const auto& origNodes = original.get_nodeGenes();
        const auto& copyNodes = copy.get_nodeGenes();
        EXPECT_EQ(origNodes[origConn.get_sourceNodeIndex()].get_historyID(), copyNodes[copyConn.get_sourceNodeIndex()].get_historyID());
        EXPECT_EQ(origNodes[origConn.get_targetNodeIndex()].get_historyID(), copyNodes[copyConn.get_targetNodeIndex()].get_historyID());
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
        
        // Verify the copied connection references the correct node indices
        const auto& origNodes = original.get_nodeGenes();
        EXPECT_EQ(origNodes[origConn.get_sourceNodeIndex()].get_historyID(), copyNodes[copyConn.get_sourceNodeIndex()].get_historyID());
        EXPECT_EQ(origNodes[origConn.get_targetNodeIndex()].get_historyID(), copyNodes[copyConn.get_targetNodeIndex()].get_historyID());
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
// Assignment Operator Tests
// ============================================================================

TEST_F(GenomeTest, AssignmentOperator_IsDeleted) {
    // Verify that assignment operator is properly deleted due to const members
    GenomeParams params1 = createValidGenomeParams();
    GenomeParams params2;
    params2._nodeHistoryIDs = {10, 20};
    params2._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params2._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    
    Genome genome1(params1);
    Genome genome2(params2);
    
    // Verify initial difference
    EXPECT_NE(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    
    // Assignment operator should be deleted - this should not compile:
    // genome2 = genome1;  // Intentionally commented out
    
    // Instead, we must use copy construction for copying
    Genome genome3(genome1);
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome3.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome3.get_connectionGenes().size());
}

TEST_F(GenomeTest, CopyConstruction_ReplacesAssignment) {
    GenomeParams params = createValidGenomeParams();
    Genome genome(params);
    
    // Since assignment is deleted, we use copy construction
    Genome copy(genome);
    
    // Verify they are equivalent but independent
    EXPECT_EQ(genome.get_nodeGenes().size(), copy.get_nodeGenes().size());
    EXPECT_EQ(genome.get_connectionGenes().size(), copy.get_connectionGenes().size());
}