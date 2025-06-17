#include <gtest/gtest.h>
#include <memory>

#include "tests/test_common.h"
#include "version3/operator/CycleDetection.hpp"
#include "version3/operator/Init.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Operator;

class CycleDetectionTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_shared<HistoryTracker>();
    }

    std::shared_ptr<HistoryTracker> historyTracker;
    
    // Helper to create minimal NEAT starting genome (I/O/BIAS, linear connection)
    Genome createMinimalGenome() {
        std::vector<NodeGeneAttributes> inputAttrs = {{ActivationType::NONE}};
        std::vector<NodeGeneAttributes> outputAttrs = {{ActivationType::SIGMOID}};
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
        
        InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                         InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
        
        return init(historyTracker, params);
    }
    
    // Helper to create genome with specific nodes and connections
    Genome createCustomGenome(const std::vector<uint32_t>& nodeHistoryIDs,
                             const std::vector<NodeType>& nodeTypes,
                             const std::vector<NodeGeneAttributes>& nodeAttrs,
                             const std::vector<uint32_t>& connHistoryIDs,
                             const std::vector<uint32_t>& sourceIDs,
                             const std::vector<uint32_t>& targetIDs,
                             const std::vector<ConnectionGeneAttributes>& connAttrs) {
        GenomeParams params;
        params._nodeHistoryIDs = nodeHistoryIDs;
        params._nodeTypes = nodeTypes;
        params._nodeAttributes = nodeAttrs;
        params._connectionHistoryIDs = connHistoryIDs;
        params._sourceNodeHistoryIDs = sourceIDs;
        params._targetNodeHistoryIDs = targetIDs;
        params._connectionAttributes = connAttrs;
        
        return Genome(params);
    }
};

// =============================================================================
// 1. CORE CYCLE DETECTION FUNCTIONALITY
// =============================================================================

TEST_F(CycleDetectionTest, EmptyGenome) {
    // Empty genome should have no cycles
    GenomeParams params;
    Genome emptyGenome(params);
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(emptyGenome, cycleParams)) << "Empty genome should have no cycles";
}

TEST_F(CycleDetectionTest, SingleIsolatedNode) {
    // Single node with no connections should have no cycles
    Genome genome = createCustomGenome(
        {1}, {NodeType::HIDDEN}, {{ActivationType::SIGMOID}},
        {}, {}, {}, {}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Single isolated node should have no cycles";
}

TEST_F(CycleDetectionTest, LinearChain_NoCycles) {
    // Linear chain A→B→C should have no cycles
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 3},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Linear chain A→B→C should have no cycles";
}

TEST_F(CycleDetectionTest, MultipleDisconnectedLinearChains) {
    // Two separate linear chains: A→B and C→D should have no cycles
    Genome genome = createCustomGenome(
        {1, 2, 3, 4}, 
        {NodeType::INPUT, NodeType::HIDDEN, NodeType::INPUT, NodeType::HIDDEN},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 3}, 
        {2, 4},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Multiple disconnected linear chains should have no cycles";
}

TEST_F(CycleDetectionTest, SelfLoop) {
    // Self-loop: A→A should detect cycle
    Genome genome = createCustomGenome(
        {1}, 
        {NodeType::HIDDEN},
        {{ActivationType::SIGMOID}},
        {1}, 
        {1}, 
        {1},
        {{1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Self-loop A→A should detect cycle";
}

TEST_F(CycleDetectionTest, TwoNodeCycle) {
    // Two-node cycle: A→B→A should detect cycle
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Two-node cycle A→B→A should detect cycle";
}

TEST_F(CycleDetectionTest, ThreeNodeCycle) {
    // Three-node cycle: A→B→C→A should detect cycle
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3}, 
        {1, 2, 3}, 
        {2, 3, 1},
        {{1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Three-node cycle A→B→C→A should detect cycle";
}

// =============================================================================
// 2. CONNECTION STATE HANDLING (CRITICAL FOR NEAT)
// =============================================================================

TEST_F(CycleDetectionTest, CycleExistsOnlyThroughDisabledConnections) {
    // Cycle A→B→A where one connection is disabled should return false
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, true}, {1.0f, false}}  // Second connection disabled
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Cycle with disabled connection should return false";
}

TEST_F(CycleDetectionTest, CycleExistsThroughEnabledConnections) {
    // Cycle A→B→A where both connections are enabled should return true
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, true}, {1.0f, true}}  // Both connections enabled
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Cycle with enabled connections should return true";
}

TEST_F(CycleDetectionTest, CycleBrokenByDisablingConnection) {
    // Three-node cycle A→B→C→A with middle connection disabled
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3}, 
        {1, 2, 3}, 
        {2, 3, 1},
        {{1.0f, true}, {1.0f, false}, {1.0f, true}}  // Middle connection disabled
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Cycle broken by disabled connection should return false";
}

// =============================================================================
// 3. NODE TYPE FILTERING PARAMETER
// =============================================================================

TEST_F(CycleDetectionTest, IncludeAllNodeTypes_CycleWithHiddenAndBias) {
    // Cycle involving HIDDEN→BIAS→HIDDEN should be detected when includeAllNodeTypes=true
    // Note: INPUT nodes cannot be targets, OUTPUT nodes cannot be sources
    Genome genome = createCustomGenome(
        {1, 2, 3, 4}, 
        {NodeType::INPUT, NodeType::HIDDEN, NodeType::BIAS, NodeType::OUTPUT},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}, {ActivationType::SIGMOID}},
        {1, 2, 3, 4}, 
        {1, 2, 3, 2}, 
        {2, 3, 2, 4},  // INPUT→HIDDEN, HIDDEN→BIAS, BIAS→HIDDEN (cycle), HIDDEN→OUTPUT
        {{1.0f, true}, {1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(true);  // includeAllNodeTypes = true
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Cycle with HIDDEN→BIAS→HIDDEN should be detected when includeAllNodeTypes=true";
}

TEST_F(CycleDetectionTest, IncludeAllNodeTypes_CycleWithBiasHidden) {
    // Cycle involving BIAS→HIDDEN→BIAS should be detected when includeAllNodeTypes=true
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::BIAS, NodeType::HIDDEN},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(true);  // includeAllNodeTypes = true
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Cycle with BIAS→HIDDEN should be detected when includeAllNodeTypes=true";
}

TEST_F(CycleDetectionTest, HiddenNodesOnly_CycleWithOnlyHidden) {
    // Cycle involving only HIDDEN nodes should be detected when includeAllNodeTypes=false
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(false);  // includeAllNodeTypes = false (only HIDDEN)
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Cycle with only HIDDEN nodes should be detected when includeAllNodeTypes=false";
}

TEST_F(CycleDetectionTest, HiddenNodesOnly_CycleWithoutHiddenNodes) {
    // Cycle involving only INPUT→BIAS (no HIDDEN nodes) should NOT be detected when includeAllNodeTypes=false
    // Note: INPUT cannot be target, OUTPUT cannot be source, so only valid flow is INPUT→BIAS/OUTPUT
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::INPUT, NodeType::BIAS, NodeType::OUTPUT},
        {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 3},  // INPUT→BIAS→OUTPUT (no cycle, no HIDDEN nodes)
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(false);  // includeAllNodeTypes = false (only HIDDEN)
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Genome with no HIDDEN nodes should NOT detect cycles when includeAllNodeTypes=false";
}

TEST_F(CycleDetectionTest, HiddenNodesOnly_HiddenBiasHiddenFiltered) {
    // HIDDEN→BIAS→HIDDEN where BIAS is filtered out should NOT detect cycle
    // Note: INPUT cannot be target, so using BIAS instead
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::HIDDEN, NodeType::BIAS, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::NONE}, {ActivationType::SIGMOID}},
        {1, 2, 3}, 
        {1, 2, 3}, 
        {2, 3, 1},  // Cycle: HIDDEN(1)→BIAS(2)→HIDDEN(3)→HIDDEN(1), broken when BIAS filtered
        {{1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(false);  // includeAllNodeTypes = false (only HIDDEN)
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "HIDDEN→BIAS→HIDDEN should NOT detect cycle when BIAS filtered out";
}

// =============================================================================
// 4. EDGE CASES THAT COULD OCCUR IN NEAT
// =============================================================================

TEST_F(CycleDetectionTest, NoNodesMatchFilteringCriteria) {
    // Genome with no HIDDEN nodes when includeAllNodeTypes=false
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}},
        {1}, 
        {1}, 
        {2},
        {{1.0f, true}}
    );
    
    CycleDetectionParams cycleParams(false);  // includeAllNodeTypes = false (only HIDDEN)
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Genome with no HIDDEN nodes should return false when includeAllNodeTypes=false";
}

TEST_F(CycleDetectionTest, AllConnectionsDisabled) {
    // Genome with potential cycle but all connections disabled
    Genome genome = createCustomGenome(
        {1, 2}, 
        {NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 2}, 
        {2, 1},
        {{1.0f, false}, {1.0f, false}}  // All connections disabled
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Genome with all connections disabled should have no cycles";
}

TEST_F(CycleDetectionTest, DisconnectedComponents_CycleInOneComponentOnly) {
    // Cycle in one component (nodes 1-2), acyclic in another (nodes 3-4)
    Genome genome = createCustomGenome(
        {1, 2, 3, 4}, 
        {NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3}, 
        {1, 2, 3}, 
        {2, 1, 4},  // Component 1: 1→2→1 (cycle), Component 2: 3→4 (acyclic)
        {{1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Should detect cycle even when only one component has cycles";
}

// =============================================================================
// 5. REALISTIC NEAT GENOME SCENARIOS
// =============================================================================

TEST_F(CycleDetectionTest, TypicalMinimalNEATStartingGenome) {
    // Typical minimal NEAT starting genome should be acyclic
    Genome genome = createMinimalGenome();
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Typical minimal NEAT starting genome should be acyclic";
}

TEST_F(CycleDetectionTest, PostNodeMutationGenome) {
    // Genome after node mutation (split connection) should remain acyclic
    // Original: INPUT→OUTPUT, After split: INPUT→HIDDEN→OUTPUT
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2}, 
        {1, 3}, 
        {3, 2},
        {{1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_FALSE(hasCycles(genome, cycleParams)) << "Genome after node mutation should remain acyclic";
}

TEST_F(CycleDetectionTest, PostConnectionMutationGenome_CreatesCycle) {
    // Genome after connection mutation that creates a cycle
    // Linear: 1→2→3, with added connection: 3→2 (creates cycle between HIDDEN nodes)
    Genome genome = createCustomGenome(
        {1, 2, 3}, 
        {NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3}, 
        {1, 2, 3}, 
        {2, 3, 2},  // Connections: 1→2, 2→3, 3→2 (creates cycle 2↔3)
        {{1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Genome after connection mutation could create cycle";
}

TEST_F(CycleDetectionTest, PostCrossoverGenome_AccidentalCycle) {
    // Crossover result that accidentally creates cycle through inherited connections
    // Parent A: 1→2, 3→1  Parent B: 2→3, 1→4  Offspring: 1→2→3→1 (cycle in HIDDEN layer)
    Genome genome = createCustomGenome(
        {1, 2, 3, 4}, 
        {NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3, 4}, 
        {1, 2, 3, 1}, 
        {2, 3, 1, 4},  // Creates cycle: 1→2→3→1, plus output 1→4
        {{1.0f, true}, {1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Crossover result could accidentally create cycle";
}

// =============================================================================
// 6. FAIL-FAST VERIFICATION (BEHAVIOR ONLY)
// =============================================================================

TEST_F(CycleDetectionTest, MultipleCyclesPresent_ReturnsTrue) {
    // Multiple cycles present: 1→2→1 and 3→4→3
    // Verify returns true immediately (don't measure timing, just verify it finds A cycle)
    Genome genome = createCustomGenome(
        {1, 2, 3, 4}, 
        {NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN},
        {{ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}, {ActivationType::SIGMOID}},
        {1, 2, 3, 4}, 
        {1, 2, 3, 4}, 
        {2, 1, 4, 3},  // Two cycles: 1→2→1 and 3→4→3
        {{1.0f, true}, {1.0f, true}, {1.0f, true}, {1.0f, true}}
    );
    
    CycleDetectionParams cycleParams;
    
    EXPECT_TRUE(hasCycles(genome, cycleParams)) << "Should return true when multiple cycles present";
}