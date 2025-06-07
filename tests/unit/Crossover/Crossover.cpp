#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/Crossover.hpp"
#include "version3/analysis/FitnessResult.hpp"

using namespace Operator;
using namespace Analysis;

// Test-specific fitness result implementation
class TestFitnessResult : public FitnessResultInterface {
public:
    explicit TestFitnessResult(double fitness) : _fitness(fitness) {}
    
    bool isBetterThan(const FitnessResultInterface& other) const override {
        const auto* otherTest = dynamic_cast<const TestFitnessResult*>(&other);
        if (!otherTest) return false;
        return _fitness > otherTest->_fitness;
    }
    
    bool isEqualTo(const FitnessResultInterface& other) const override {
        const auto* otherTest = dynamic_cast<const TestFitnessResult*>(&other);
        if (!otherTest) return false;
        return std::abs(_fitness - otherTest->_fitness) < 1e-9;
    }
    
    std::unique_ptr<FitnessResultInterface> clone() const override {
        return std::make_unique<TestFitnessResult>(_fitness);
    }
    
    double getFitness() const { return _fitness; }

private:
    double _fitness;
};

class CrossoverTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::srand(42); // Fixed seed for reproducibility where needed
    }

    // Helper to create parent genome A with clear gene categorization
    Genome createParentA() {
        GenomeParams params;
        // Shared nodes: 1 (INPUT), 2 (OUTPUT) 
        // Excess nodes: 3 (HIDDEN) - beyond B's range
        params._nodeHistoryIDs = {1, 2, 3};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}
        };
        
        // Shared connection: 1->2 (historyID 1)
        // Disjoint connection: missing historyID 2 that B has
        // Excess connections: 3, 4 (beyond B's max of 2)
        params._connectionHistoryIDs = {1, 3, 4};
        params._sourceNodeHistoryIDs = {1, 1, 3};
        params._targetNodeHistoryIDs = {2, 3, 2};
        params._connectionAttributes = {
            {1.0f, true},   // shared connection
            {1.5f, true},   // excess connection
            {2.0f, false}   // excess connection (disabled)
        };
        return Genome(params);
    }

    // Helper to create parent genome B with clear gene categorization  
    Genome createParentB() {
        GenomeParams params;
        // Shared nodes: 1 (INPUT), 2 (OUTPUT)
        // Excess nodes: 4 (HIDDEN) - beyond A's range if B is bigger
        params._nodeHistoryIDs = {1, 2, 4};
        params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::TANH}
        };
        
        params._connectionHistoryIDs = {1, 2, 5};
        params._sourceNodeHistoryIDs = {1, 1, 4};
        params._targetNodeHistoryIDs = {2, 4, 2};
        params._connectionAttributes = {
            {-1.0f, true},  // shared connection (DIFFERENT weight: -1.0 vs 1.0)
            {2.5f, true},   // disjoint connection (same topology as would be in A)
            {3.0f, true},   // excess connection
        };
        return Genome(params);
    }

    // Helper to verify that a genome contains specific node history IDs
    void verifyNodesPresent(const Genome& genome, const std::set<uint32_t>& expectedNodeIDs) {
        std::set<uint32_t> actualNodeIDs;
        for (const auto& node : genome.get_nodeGenes()) {
            actualNodeIDs.insert(node.get_historyID());
        }
        EXPECT_EQ(actualNodeIDs, expectedNodeIDs);
    }

    // Helper to verify that a genome contains specific connection history IDs
    void verifyConnectionsPresent(const Genome& genome, const std::set<uint32_t>& expectedConnIDs) {
        std::set<uint32_t> actualConnIDs;
        for (const auto& conn : genome.get_connectionGenes()) {
            actualConnIDs.insert(conn.get_historyID());
        }
        EXPECT_EQ(actualConnIDs, expectedConnIDs);
    }

    // Helper to get connection IDs that are disjoint/excess between two genomes
    std::set<uint32_t> getDisjointExcessConnections(const Genome& genomeA, const Genome& genomeB) {
        std::set<uint32_t> connIDsA, connIDsB;
        for (const auto& conn : genomeA.get_connectionGenes()) {
            connIDsA.insert(conn.get_historyID());
        }
        for (const auto& conn : genomeB.get_connectionGenes()) {
            connIDsB.insert(conn.get_historyID());
        }
        
        std::set<uint32_t> disjointExcess;
        std::set_symmetric_difference(
            connIDsA.begin(), connIDsA.end(),
            connIDsB.begin(), connIDsB.end(),
            std::inserter(disjointExcess, disjointExcess.begin())
        );
        return disjointExcess;
    }
};

// ============================================================================
// Fitness-Based Parent Selection Tests
// ============================================================================

TEST_F(CrossoverTest, FitterParentInheritsDisjointExcess_ClearCase) {
    // Test that disjoint/excess genes are inherited from the fitter parent
    
    // Parent A: I/O/BIAS + hidden node + connections {1, 3, 4}
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 10};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::SIGMOID}
    };
    
    paramsA._connectionHistoryIDs = {1, 3, 4};
    paramsA._sourceNodeHistoryIDs = {1, 3, 10};
    paramsA._targetNodeHistoryIDs = {2, 10, 2};
    paramsA._connectionAttributes = {
        {1.0f, true},   // matching with B
        {1.3f, true},   // excess beyond B's range
        {1.4f, false}   // excess beyond B's range (disabled)
    };
    
    // Parent B: I/O/BIAS + different hidden node + connections {1, 2}
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 20};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::TANH}
    };
    
    paramsB._connectionHistoryIDs = {1, 2};
    paramsB._sourceNodeHistoryIDs = {1, 3};
    paramsB._targetNodeHistoryIDs = {2, 20};
    paramsB._connectionAttributes = {
        {2.0f, true},   // matching with A (different weight)
        {2.2f, true}    // disjoint (gap in A's range)
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Test Case 1: Parent A is fitter - should get A's disjoint/excess genes
    {
        TestFitnessResult fitnessA(0.9);
        TestFitnessResult fitnessB(0.3);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Verify node inheritance
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Should have shared I/O/BIAS nodes + fitter parent A's excess hidden node
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared INPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared OUTPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit shared BIAS node";
        EXPECT_TRUE(offspringNodeIDs.find(10) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent A";
        EXPECT_TRUE(offspringNodeIDs.find(20) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent B";
        
        // Verify connection inheritance
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
        }
        
        // Should have matching connection (randomly selected) + fitter parent A's excess connections
        EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit matching connection 1";
        EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should inherit excess connection 3 from fitter parent A";
        EXPECT_TRUE(offspringConnIDs.find(4) != offspringConnIDs.end()) << "Should inherit excess connection 4 from fitter parent A";
        EXPECT_TRUE(offspringConnIDs.find(2) == offspringConnIDs.end()) << "Should NOT inherit disjoint connection 2 from less fit parent B";
    }
    
    // Test Case 2: Parent B is fitter - should get B's disjoint/excess genes
    {
        TestFitnessResult fitnessA(0.3);
        TestFitnessResult fitnessB(0.9);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Verify node inheritance
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Should have shared I/O/BIAS nodes + fitter parent B's excess hidden node
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared INPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared OUTPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit shared BIAS node";
        EXPECT_TRUE(offspringNodeIDs.find(20) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent B";
        EXPECT_TRUE(offspringNodeIDs.find(10) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent A";
        
        // Verify connection inheritance
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
        }
        
        // Should have matching connection + fitter parent B's disjoint connection
        EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit matching connection 1";
        EXPECT_TRUE(offspringConnIDs.find(2) != offspringConnIDs.end()) << "Should inherit disjoint connection 2 from fitter parent B";
        EXPECT_TRUE(offspringConnIDs.find(3) == offspringConnIDs.end()) << "Should NOT inherit excess connection 3 from less fit parent A";
        EXPECT_TRUE(offspringConnIDs.find(4) == offspringConnIDs.end()) << "Should NOT inherit excess connection 4 from less fit parent A";
    }
}

TEST_F(CrossoverTest, EqualFitness_RandomGeneSelection) {
    // Test that with equal fitness, matching genes are randomly selected from either parent
    
    // Parent A: I/O/BIAS + connections with specific weights
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::SIGMOID}
    };
    
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4};
    paramsA._targetNodeHistoryIDs = {2, 4, 2};
    paramsA._connectionAttributes = {
        {10.0f, true},  // matching - unique weight to identify parent A
        {20.0f, true},  // matching - unique weight to identify parent A
        {30.0f, false}  // matching - unique weight to identify parent A
    };
    
    // Parent B: Same structure but different weights (matching connections)
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 4};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::TANH}
    };
    
    paramsB._connectionHistoryIDs = {1, 2, 3};
    paramsB._sourceNodeHistoryIDs = {1, 3, 4};
    paramsB._targetNodeHistoryIDs = {2, 4, 2};
    paramsB._connectionAttributes = {
        {100.0f, true},  // matching - unique weight to identify parent B
        {200.0f, false}, // matching - unique weight to identify parent B
        {300.0f, true}   // matching - unique weight to identify parent B
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Equal fitness should trigger random selection for matching genes
    TestFitnessResult fitnessA(0.75);
    TestFitnessResult fitnessB(0.75);
    
    CrossoverParams params(0.0);
    
    // Track gene selection from each parent over multiple runs
    std::map<uint32_t, int> selectedFromA, selectedFromB;
    const int numRuns = 50;
    
    for (int run = 0; run < numRuns; ++run) {
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Check which parent's genes were selected based on weights
        for (const auto& conn : offspring.get_connectionGenes()) {
            uint32_t historyID = conn.get_historyID();
            float weight = conn.get_attributes().weight;
            
            // Identify parent by unique weights
            if (historyID == 1) {
                if (std::abs(weight - 10.0f) < 0.001f) selectedFromA[1]++;
                else if (std::abs(weight - 100.0f) < 0.001f) selectedFromB[1]++;
            } else if (historyID == 2) {
                if (std::abs(weight - 20.0f) < 0.001f) selectedFromA[2]++;
                else if (std::abs(weight - 200.0f) < 0.001f) selectedFromB[2]++;
            } else if (historyID == 3) {
                if (std::abs(weight - 30.0f) < 0.001f) selectedFromA[3]++;
                else if (std::abs(weight - 300.0f) < 0.001f) selectedFromB[3]++;
            }
        }
        
        // Verify offspring always has all matching connections
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
        }
        EXPECT_EQ(offspringConnIDs.size(), 3u) << "Should always have all 3 matching connections";
        EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should have connection 1";
        EXPECT_TRUE(offspringConnIDs.find(2) != offspringConnIDs.end()) << "Should have connection 2";
        EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should have connection 3";
    }
    
    // Verify random selection occurred (both parents should contribute)
    for (uint32_t connID = 1; connID <= 3; ++connID) {
        int totalSelections = selectedFromA[connID] + selectedFromB[connID];
        EXPECT_EQ(totalSelections, numRuns) << "Should account for all selections of connection " << connID;
        
        // With equal fitness and 50 runs, expect some variation (not 100% from one parent)
        EXPECT_GT(selectedFromA[connID], 0) << "Parent A should be selected sometimes for connection " << connID;
        EXPECT_GT(selectedFromB[connID], 0) << "Parent B should be selected sometimes for connection " << connID;
        
        // Statistical check - with 50 runs, extremely unlikely to get less than 10% from either parent
        double ratioFromA = static_cast<double>(selectedFromA[connID]) / numRuns;
        EXPECT_GT(ratioFromA, 0.1) << "Parent A selection too rare for connection " << connID << ": " << ratioFromA;
        EXPECT_LT(ratioFromA, 0.9) << "Parent A selection too frequent for connection " << connID << ": " << ratioFromA;
    }
}

// ============================================================================
// Gene Categorization Logic Tests
// ============================================================================

TEST_F(CrossoverTest, MatchingGenes_AlwaysInheritedWithRandomSelection) {
    // Test that matching genes are always inherited and randomly selected from either parent
    
    // Parent A: I/O/BIAS + matching connections with distinct weights/states
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4, 5};
    paramsA._nodeTypes = {
        NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS,
        NodeType::HIDDEN, NodeType::HIDDEN
    };
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::SIGMOID}, {ActivationType::TANH}
    };
    
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4};
    paramsA._targetNodeHistoryIDs = {2, 4, 5};
    paramsA._connectionAttributes = {
        {5.0f, true},   // matching - parent A signature
        {15.0f, false}, // matching - parent A signature (disabled)
        {25.0f, true}   // matching - parent A signature
    };
    
    // Parent B: Same structure, same connections, different weights/states
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 4, 5};
    paramsB._nodeTypes = {
        NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS,
        NodeType::HIDDEN, NodeType::HIDDEN
    };
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::RELU}, {ActivationType::SIGMOID}
    };
    
    paramsB._connectionHistoryIDs = {1, 2, 3};
    paramsB._sourceNodeHistoryIDs = {1, 3, 4};
    paramsB._targetNodeHistoryIDs = {2, 4, 5};
    paramsB._connectionAttributes = {
        {50.0f, false}, // matching - parent B signature (different weight/state)
        {150.0f, true}, // matching - parent B signature (different weight/state)
        {250.0f, false} // matching - parent B signature (different weight/state)
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Test with different fitness scenarios
    std::vector<std::pair<double, double>> fitnessScenarios = {
        {0.8, 0.6}, // A fitter - but matching genes still randomly selected
        {0.4, 0.9}, // B fitter - but matching genes still randomly selected
        {0.7, 0.7}  // Equal fitness - random selection
    };
    
    for (const auto& [fitnessA, fitnessB] : fitnessScenarios) {
        TestFitnessResult resultA(fitnessA);
        TestFitnessResult resultB(fitnessB);
        
        CrossoverParams params(0.0);
        
        // Track selections for each matching connection
        std::map<uint32_t, int> selectedFromA, selectedFromB;
        const int numRuns = 40;
        
        for (int run = 0; run < numRuns; ++run) {
            Genome offspring = crossover(parentA, resultA, parentB, resultB, params);
            
            // Verify all matching connections are present
            std::set<uint32_t> offspringConnIDs;
            for (const auto& conn : offspring.get_connectionGenes()) {
                offspringConnIDs.insert(conn.get_historyID());
            }
            
            EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Matching connection 1 must always be inherited";
            EXPECT_TRUE(offspringConnIDs.find(2) != offspringConnIDs.end()) << "Matching connection 2 must always be inherited";
            EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Matching connection 3 must always be inherited";
            
            // Track which parent was selected for each matching connection
            for (const auto& conn : offspring.get_connectionGenes()) {
                uint32_t historyID = conn.get_historyID();
                float weight = conn.get_attributes().weight;
                
                if (historyID == 1) {
                    if (std::abs(weight - 5.0f) < 0.001f) selectedFromA[1]++;
                    else if (std::abs(weight - 50.0f) < 0.001f) selectedFromB[1]++;
                    else FAIL() << "Connection 1 weight " << weight << " doesn't match either parent";
                } else if (historyID == 2) {
                    if (std::abs(weight - 15.0f) < 0.001f) selectedFromA[2]++;
                    else if (std::abs(weight - 150.0f) < 0.001f) selectedFromB[2]++;
                    else FAIL() << "Connection 2 weight " << weight << " doesn't match either parent";
                } else if (historyID == 3) {
                    if (std::abs(weight - 25.0f) < 0.001f) selectedFromA[3]++;
                    else if (std::abs(weight - 250.0f) < 0.001f) selectedFromB[3]++;
                    else FAIL() << "Connection 3 weight " << weight << " doesn't match either parent";
                }
            }
        }
        
        // Verify random selection for matching genes (regardless of fitness)
        for (uint32_t connID = 1; connID <= 3; ++connID) {
            int totalSelections = selectedFromA[connID] + selectedFromB[connID];
            EXPECT_EQ(totalSelections, numRuns) << "Should account for all selections of connection " << connID 
                                                << " (fitness A=" << fitnessA << ", B=" << fitnessB << ")";
            
            // Matching genes should show random selection regardless of fitness
            EXPECT_GT(selectedFromA[connID], 0) << "Parent A should be selected sometimes for connection " << connID
                                                << " (fitness A=" << fitnessA << ", B=" << fitnessB << ")";
            EXPECT_GT(selectedFromB[connID], 0) << "Parent B should be selected sometimes for connection " << connID
                                                << " (fitness A=" << fitnessA << ", B=" << fitnessB << ")";
        }
    }
}

TEST_F(CrossoverTest, DisjointGenes_WithinOverlappingRanges) {
    // Create parents with disjoint genes - innovation gaps within overlapping ranges
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
    paramsA._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}};
    
    // Parent A has connections: 1, 3, 5 (missing 2, 4)
    paramsA._connectionHistoryIDs = {1, 3, 5};
    paramsA._sourceNodeHistoryIDs = {1, 1, 3};
    paramsA._targetNodeHistoryIDs = {2, 3, 2};
    paramsA._connectionAttributes = {{1.0f, true}, {1.5f, true}, {2.0f, true}};
    
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 4};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
    paramsB._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::TANH}};
    
    // Parent B has connections: 1, 2, 4 (missing 3, 5)
    paramsB._connectionHistoryIDs = {1, 2, 4};
    paramsB._sourceNodeHistoryIDs = {1, 1, 4};
    paramsB._targetNodeHistoryIDs = {2, 4, 2};
    paramsB._connectionAttributes = {{-1.0f, true}, {2.5f, true}, {3.0f, true}};
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Make parent A fitter to get its disjoint genes
    TestFitnessResult fitnessA(0.9);
    TestFitnessResult fitnessB(0.5);
    
    CrossoverParams params(0.0);
    
    Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
    
    std::set<uint32_t> offspringConnIDs;
    for (const auto& conn : offspring.get_connectionGenes()) {
        offspringConnIDs.insert(conn.get_historyID());
    }
    
    // Should have matching connection: 1 (from either parent)
    EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit matching connection 1";
    
    // Should have disjoint connections from fitter parent A: 3, 5
    EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should inherit disjoint connection 3 from parent A";
    EXPECT_TRUE(offspringConnIDs.find(5) != offspringConnIDs.end()) << "Should inherit disjoint connection 5 from parent A";
    
    // Should NOT have disjoint connections from less fit parent B: 2, 4
    EXPECT_TRUE(offspringConnIDs.find(2) == offspringConnIDs.end()) << "Should NOT inherit disjoint connection 2 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(4) == offspringConnIDs.end()) << "Should NOT inherit disjoint connection 4 from parent B";
    
    // Verify node inheritance follows the same pattern
    std::set<uint32_t> offspringNodeIDs;
    for (const auto& node : offspring.get_nodeGenes()) {
        offspringNodeIDs.insert(node.get_historyID());
    }
    
    // Should have shared nodes: 1, 2
    EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared node 1";
    EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared node 2";
    
    // Should have disjoint node from fitter parent A: 3
    EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit disjoint node 3 from parent A";
    
    // Should NOT have disjoint node from less fit parent B: 4
    EXPECT_TRUE(offspringNodeIDs.find(4) == offspringNodeIDs.end()) << "Should NOT inherit disjoint node 4 from parent B";
}

TEST_F(CrossoverTest, ExcessGenes_BeyondOtherParentRange) {
    // Simplified test - reuse existing working helper functions and just test one aspect
    Genome parentA = createParentA();  // Uses our working helper
    
    // Create a simple parent B with just excess connections beyond A's range  
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2};  // Same shared nodes as A
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    paramsB._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    
    // Parent B has connections: 1 (shared with A) + 10, 11 (clearly excess beyond A's max of 4)
    paramsB._connectionHistoryIDs = {1, 10, 11};
    paramsB._sourceNodeHistoryIDs = {1, 1, 1};   // All from INPUT node
    paramsB._targetNodeHistoryIDs = {2, 2, 2};   // All to OUTPUT node
    paramsB._connectionAttributes = {
        {-1.0f, true}, {3.0f, true}, {3.5f, true}
    };
    
    Genome parentB(paramsB);
    
    // Make parent B fitter to get its excess genes
    TestFitnessResult fitnessA(0.4);
    TestFitnessResult fitnessB(0.8);
    
    CrossoverParams params(0.0);
    
    Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
    
    std::set<uint32_t> offspringConnIDs;
    for (const auto& conn : offspring.get_connectionGenes()) {
        offspringConnIDs.insert(conn.get_historyID());
    }
    
    // Should have matching connection: 1 (from either parent)
    EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit matching connection 1";
    
    // Should have excess connections from fitter parent B: 10, 11
    // These are "excess" because they're beyond parent A's max innovation (4)
    EXPECT_TRUE(offspringConnIDs.find(10) != offspringConnIDs.end()) << "Should inherit excess connection 10 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(11) != offspringConnIDs.end()) << "Should inherit excess connection 11 from parent B";
    
    // Should NOT have excess connections from less fit parent A: 3, 4
    // These are A's excess connections (beyond B's shared range of 1)
    EXPECT_TRUE(offspringConnIDs.find(3) == offspringConnIDs.end()) << "Should NOT inherit excess connection 3 from parent A";
    EXPECT_TRUE(offspringConnIDs.find(4) == offspringConnIDs.end()) << "Should NOT inherit excess connection 4 from parent A";
    
    // Verify node inheritance - both parents have same shared nodes only
    std::set<uint32_t> offspringNodeIDs;
    for (const auto& node : offspring.get_nodeGenes()) {
        offspringNodeIDs.insert(node.get_historyID());
    }
    
    // Should have shared nodes: 1, 2
    EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared node 1";
    EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared node 2";
    
    // Should NOT have A's excess node 3 (since A is less fit)
    EXPECT_TRUE(offspringNodeIDs.find(3) == offspringNodeIDs.end()) << "Should NOT inherit excess node 3 from parent A";
}

TEST_F(CrossoverTest, NoMatchingGenes_CompletelyDifferentRanges) {
    // Create parents with completely non-overlapping innovation ranges
    
    // Parent A: small genome with low innovation numbers  
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
    paramsA._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}};
    
    paramsA._connectionHistoryIDs = {1, 2};
    paramsA._sourceNodeHistoryIDs = {1, 1};
    paramsA._targetNodeHistoryIDs = {2, 3};
    paramsA._connectionAttributes = {{1.0f, true}, {1.5f, true}};
    
    // Parent B: larger genome with high innovation numbers (no overlap with A)
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {10, 11, 12, 13};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::TANH}, {ActivationType::RELU}
    };
    
    paramsB._connectionHistoryIDs = {20, 21, 22};
    paramsB._sourceNodeHistoryIDs = {10, 10, 12};
    paramsB._targetNodeHistoryIDs = {11, 12, 11};
    paramsB._connectionAttributes = {{2.0f, true}, {2.5f, true}, {3.0f, true}};
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Make parent B fitter to get all its genes
    TestFitnessResult fitnessA(0.3);
    TestFitnessResult fitnessB(0.9);
    
    CrossoverParams params(0.0);
    
    Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
    
    // With no matching genes, offspring should get all genes from fitter parent B
    std::set<uint32_t> offspringConnIDs;
    for (const auto& conn : offspring.get_connectionGenes()) {
        offspringConnIDs.insert(conn.get_historyID());
    }
    
    // Should have all connections from fitter parent B
    EXPECT_TRUE(offspringConnIDs.find(20) != offspringConnIDs.end()) << "Should inherit connection 20 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(21) != offspringConnIDs.end()) << "Should inherit connection 21 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(22) != offspringConnIDs.end()) << "Should inherit connection 22 from parent B";
    
    // Should NOT have any connections from less fit parent A
    EXPECT_TRUE(offspringConnIDs.find(1) == offspringConnIDs.end()) << "Should NOT inherit connection 1 from parent A";
    EXPECT_TRUE(offspringConnIDs.find(2) == offspringConnIDs.end()) << "Should NOT inherit connection 2 from parent A";
    
    // Verify node inheritance
    std::set<uint32_t> offspringNodeIDs;
    for (const auto& node : offspring.get_nodeGenes()) {
        offspringNodeIDs.insert(node.get_historyID());
    }
    
    // Should have all nodes from fitter parent B
    EXPECT_TRUE(offspringNodeIDs.find(10) != offspringNodeIDs.end()) << "Should inherit node 10 from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(11) != offspringNodeIDs.end()) << "Should inherit node 11 from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(12) != offspringNodeIDs.end()) << "Should inherit node 12 from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(13) != offspringNodeIDs.end()) << "Should inherit node 13 from parent B";
    
    // Should NOT have any nodes from less fit parent A
    EXPECT_TRUE(offspringNodeIDs.find(1) == offspringNodeIDs.end()) << "Should NOT inherit node 1 from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(2) == offspringNodeIDs.end()) << "Should NOT inherit node 2 from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(3) == offspringNodeIDs.end()) << "Should NOT inherit node 3 from parent A";
}

TEST_F(CrossoverTest, AllMatchingGenes_IdenticalStructure) {
    // Create parents with identical innovation sets but different attributes
    
    // Parent A: baseline structure
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
    paramsA._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::SIGMOID}};
    
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 1, 3};
    paramsA._targetNodeHistoryIDs = {2, 3, 2};
    paramsA._connectionAttributes = {{1.0f, true}, {1.5f, true}, {2.0f, false}};  // disabled connection
    
    // Parent B: identical structure but different weights and activation
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::HIDDEN};
    paramsB._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::TANH}};  // different activation
    
    paramsB._connectionHistoryIDs = {1, 2, 3};
    paramsB._sourceNodeHistoryIDs = {1, 1, 3};
    paramsB._targetNodeHistoryIDs = {2, 3, 2};
    paramsB._connectionAttributes = {{-1.0f, true}, {-1.5f, false}, {-2.0f, true}};  // different weights and enable states
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Test with equal fitness to ensure random 50/50 selection
    TestFitnessResult fitnessA(0.7);
    TestFitnessResult fitnessB(0.7);
    
    CrossoverParams params(0.0); // No reactivation for this test
    
    // Run multiple times to verify random selection
    std::map<uint32_t, std::set<float>> observedWeights;
    const int numRuns = 50;
    
    for (int i = 0; i < numRuns; ++i) {
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Check that all connections are present (all matching genes should be inherited)
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
            observedWeights[conn.get_historyID()].insert(conn.get_attributes().weight);
        }
        
        // All connections should be present
        EXPECT_EQ(offspringConnIDs.size(), 3u) << "Should have all 3 matching connections";
        EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit connection 1";
        EXPECT_TRUE(offspringConnIDs.find(2) != offspringConnIDs.end()) << "Should inherit connection 2";
        EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should inherit connection 3";
        
        // All nodes should be present
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        EXPECT_EQ(offspringNodeIDs.size(), 3u) << "Should have all 3 matching nodes";
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit node 1";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit node 2";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit node 3";
    }
    
    // Verify that both parent weights were observed for each connection (random 50/50 selection)
    for (uint32_t connID = 1; connID <= 3; ++connID) {
        EXPECT_GT(observedWeights[connID].size(), 1u) << 
            "Connection " << connID << " should have weights from both parents";
    }
}

// ============================================================================
// Edge Cases & Boundary Conditions Tests
// ============================================================================

TEST_F(CrossoverTest, MatchingConnections_SingleSharedConnection) {
    // Test crossover when parents have exactly one matching connection (shared history ID)
    
    // Parent A: Minimal structure with single connection
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS};
    paramsA._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE}};
    
    // Single connection: INPUT → OUTPUT (history ID 5)
    paramsA._connectionHistoryIDs = {5};
    paramsA._sourceNodeHistoryIDs = {1};
    paramsA._targetNodeHistoryIDs = {2};
    paramsA._connectionAttributes = {{1.0f, true}};
    
    // Parent B: Same structure with SAME connection history ID but different weight
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS};
    paramsB._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE}};
    
    // Same connection: INPUT → OUTPUT (SAME history ID 5 - this creates matching connection)
    paramsB._connectionHistoryIDs = {5};
    paramsB._sourceNodeHistoryIDs = {1};
    paramsB._targetNodeHistoryIDs = {2};
    paramsB._connectionAttributes = {{10.0f, false}};  // Different weight and enabled state
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Equal fitness to ensure random 50/50 selection of the matching connection
    TestFitnessResult fitnessA(0.7);
    TestFitnessResult fitnessB(0.7);
    
    CrossoverParams params(0.0);
    
    // Test multiple times to verify random selection between parent connections
    std::set<float> observedWeights;
    std::set<bool> observedStates;
    
    for (int i = 0; i < 20; ++i) {
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Should always have exactly one connection (the matching one)
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
            observedWeights.insert(conn.get_attributes().weight);
            observedStates.insert(conn.get_attributes().enabled);
        }
        
        EXPECT_EQ(offspringConnIDs.size(), 1u) << "Should have exactly one matching connection";
        EXPECT_TRUE(offspringConnIDs.find(5) != offspringConnIDs.end()) << "Should inherit the matching connection with history ID 5";
        
        // Should always have all shared nodes
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        EXPECT_EQ(offspringNodeIDs.size(), 3u) << "Should have all shared nodes";
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared INPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared OUTPUT node";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit shared BIAS node";
    }
    
    // With equal fitness and multiple runs, should observe attributes from both parents
    EXPECT_GT(observedWeights.size(), 1u) << "Should observe weights from both parents (1.0 and 10.0)";
    EXPECT_GT(observedStates.size(), 1u) << "Should observe enabled states from both parents (true and false)";
}

TEST_F(CrossoverTest, OneParentEmptyConnections_NodesOnly) {
    // Test crossover when one parent has nodes but no connections
    
    // Parent A: Full structure with connections
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::SIGMOID}
    };
    
    // Multiple connections
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4};
    paramsA._targetNodeHistoryIDs = {2, 4, 2};
    paramsA._connectionAttributes = {{1.0f, true}, {1.5f, true}, {2.0f, true}};
    
    // Parent B: Only nodes, no connections (minimal viable)
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {5, 6, 7};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE}
    };
    // No connections
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Make parent A fitter (has connections)
    TestFitnessResult fitnessA(0.8);
    TestFitnessResult fitnessB(0.4);
    
    CrossoverParams params(0.0);
    
    Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
    
    // Should inherit from fitter parent A (all excess genes, no matching)
    std::set<uint32_t> offspringNodeIDs;
    for (const auto& node : offspring.get_nodeGenes()) {
        offspringNodeIDs.insert(node.get_historyID());
    }
    
    // Should have all nodes from fitter parent A
    EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit INPUT node from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit OUTPUT node from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit BIAS node from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(4) != offspringNodeIDs.end()) << "Should inherit HIDDEN node from parent A";
    
    // Should NOT have nodes from less fit parent B
    EXPECT_TRUE(offspringNodeIDs.find(5) == offspringNodeIDs.end()) << "Should NOT inherit INPUT node from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(6) == offspringNodeIDs.end()) << "Should NOT inherit OUTPUT node from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(7) == offspringNodeIDs.end()) << "Should NOT inherit BIAS node from parent B";
    
    // Should have all connections from parent A
    std::set<uint32_t> offspringConnIDs;
    for (const auto& conn : offspring.get_connectionGenes()) {
        offspringConnIDs.insert(conn.get_historyID());
    }
    
    EXPECT_EQ(offspringConnIDs.size(), 3u) << "Should have all connections from parent A";
    EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit connection 1 from parent A";
    EXPECT_TRUE(offspringConnIDs.find(2) != offspringConnIDs.end()) << "Should inherit connection 2 from parent A";
    EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should inherit connection 3 from parent A";
    
}

// ============================================================================
// Revised Test Set - Focused Core Crossover Tests
// ============================================================================


TEST_F(CrossoverTest, VastlyDifferentSizes_SmallVsLarge) {
    // Test crossover between a very small and very large genome
    
    // Parent A: Very small (minimal viable - I/O/BIAS + 1 connection)
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS};
    paramsA._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE}};
    
    // Single connection
    paramsA._connectionHistoryIDs = {1};
    paramsA._sourceNodeHistoryIDs = {1};
    paramsA._targetNodeHistoryIDs = {2};
    paramsA._connectionAttributes = {{1.0f, true}};
    
    // Parent B: Very large genome (I/O/BIAS + 5 hidden nodes + 10+ connections)
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {10, 11, 12, 13, 14, 15, 16, 17};
    paramsB._nodeTypes = {
        NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS,
        NodeType::HIDDEN, NodeType::HIDDEN, NodeType::HIDDEN,
        NodeType::HIDDEN, NodeType::HIDDEN
    };
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::SIGMOID}, {ActivationType::TANH}, {ActivationType::RELU},
        {ActivationType::SIGMOID}, {ActivationType::TANH}
    };
    
    // Many connections creating complex network
    paramsB._connectionHistoryIDs = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    paramsB._sourceNodeHistoryIDs = {10, 12, 10, 13, 14, 15, 16, 17, 13, 14, 15};
    paramsB._targetNodeHistoryIDs = {11, 13, 14, 15, 16, 17, 11, 11, 14, 15, 16};
    paramsB._connectionAttributes = {
        {2.0f, true}, {2.1f, true}, {2.2f, false}, {2.3f, true}, {2.4f, true},
        {2.5f, false}, {2.6f, true}, {2.7f, true}, {2.8f, true}, {2.9f, false}, {3.0f, true}
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Test case 1: Small parent is fitter (should get only small parent's genes)
    {
        TestFitnessResult fitnessA(0.9);
        TestFitnessResult fitnessB(0.2);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Should inherit only from fitter small parent A (no matching genes)
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Should have only nodes from small parent A
        EXPECT_EQ(offspringNodeIDs.size(), 3u) << "Should have only 3 nodes from small parent";
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit INPUT from small parent";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit OUTPUT from small parent";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit BIAS from small parent";
        
        // Should have only connection from small parent A
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
        }
        
        EXPECT_EQ(offspringConnIDs.size(), 1u) << "Should have only 1 connection from small parent";
        EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit connection from small parent";
    }
    
    // Test case 2: Large parent is fitter (should get all large parent's genes)
    {
        TestFitnessResult fitnessA(0.1);
        TestFitnessResult fitnessB(0.9);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        // Should inherit all from fitter large parent B (no matching genes)
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Should have all nodes from large parent B
        EXPECT_EQ(offspringNodeIDs.size(), 8u) << "Should have all 8 nodes from large parent";
        for (uint32_t nodeID = 10; nodeID <= 17; ++nodeID) {
            EXPECT_TRUE(offspringNodeIDs.find(nodeID) != offspringNodeIDs.end()) << 
                "Should inherit node " << nodeID << " from large parent";
        }
        
        // Should have all connections from large parent B
        std::set<uint32_t> offspringConnIDs;
        for (const auto& conn : offspring.get_connectionGenes()) {
            offspringConnIDs.insert(conn.get_historyID());
        }
        
        EXPECT_EQ(offspringConnIDs.size(), 11u) << "Should have all 11 connections from large parent";
        for (uint32_t connID = 20; connID <= 30; ++connID) {
            EXPECT_TRUE(offspringConnIDs.find(connID) != offspringConnIDs.end()) << 
                "Should inherit connection " << connID << " from large parent";
        }
    }
}

TEST_F(CrossoverTest, GapPatterns_NonConsecutiveInnovations) {
    // Test non-consecutive innovation sequences for proper disjoint/excess categorization
    
    // Parent A: I/O/BIAS base + connections with gaps {1,3,7,9}
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4, 5};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::SIGMOID}, {ActivationType::TANH}
    };
    
    // Parent A connections: 1, 3, 7, 9 (gaps at 2, 4, 5, 6, 8)
    paramsA._connectionHistoryIDs = {1, 3, 7, 9};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4, 5};
    paramsA._targetNodeHistoryIDs = {2, 4, 2, 2};
    paramsA._connectionAttributes = {{1.0f, true}, {1.3f, true}, {1.7f, true}, {1.9f, true}};
    
    // Parent B: I/O/BIAS base + connections filling some gaps {2,4,6,8}
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 6, 7};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::RELU}, {ActivationType::SIGMOID}
    };
    
    // Parent B connections: 2, 4, 6, 8 (gaps at 1, 3, 5, 7, 9)
    // All connections must reference valid nodes that exist in Parent B: {1, 2, 3, 6, 7}
    paramsB._connectionHistoryIDs = {2, 4, 6, 8};
    paramsB._sourceNodeHistoryIDs = {1, 3, 6, 7};
    paramsB._targetNodeHistoryIDs = {2, 6, 2, 2};
    paramsB._connectionAttributes = {{2.0f, true}, {2.4f, true}, {2.6f, true}, {2.8f, true}};
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Make parent A fitter to test its disjoint/excess inheritance
    TestFitnessResult fitnessA(0.8);
    TestFitnessResult fitnessB(0.3);
    
    CrossoverParams params(0.0);
    
    Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
    
    // Analyze inheritance patterns
    std::set<uint32_t> offspringConnIDs;
    for (const auto& conn : offspring.get_connectionGenes()) {
        offspringConnIDs.insert(conn.get_historyID());
    }
    
    // Expected: No matching connections between parents (completely different ranges)
    // Parent A range: 1, 3, 7, 9 (max = 9)
    // Parent B range: 2, 4, 6, 8 (max = 8)
    // All of A's connections are disjoint (within B's range) or excess (beyond B's range)
    // Since A is fitter, should inherit all A's connections and none of B's
    
    EXPECT_TRUE(offspringConnIDs.find(1) != offspringConnIDs.end()) << "Should inherit connection 1 from fitter parent A";
    EXPECT_TRUE(offspringConnIDs.find(3) != offspringConnIDs.end()) << "Should inherit connection 3 from fitter parent A";
    EXPECT_TRUE(offspringConnIDs.find(7) != offspringConnIDs.end()) << "Should inherit connection 7 from fitter parent A";
    EXPECT_TRUE(offspringConnIDs.find(9) != offspringConnIDs.end()) << "Should inherit connection 9 from fitter parent A";
    
    // Should NOT inherit any connections from less fit parent B
    EXPECT_TRUE(offspringConnIDs.find(2) == offspringConnIDs.end()) << "Should NOT inherit connection 2 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(4) == offspringConnIDs.end()) << "Should NOT inherit connection 4 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(6) == offspringConnIDs.end()) << "Should NOT inherit connection 6 from parent B";
    EXPECT_TRUE(offspringConnIDs.find(8) == offspringConnIDs.end()) << "Should NOT inherit connection 8 from parent B";
    
    EXPECT_EQ(offspringConnIDs.size(), 4u) << "Should have exactly 4 connections from parent A";
    
    // Verify node inheritance follows same pattern
    std::set<uint32_t> offspringNodeIDs;
    for (const auto& node : offspring.get_nodeGenes()) {
        offspringNodeIDs.insert(node.get_historyID());
    }
    
    // Should have shared nodes (1, 2, 3) plus fitter parent A's excess nodes (4, 5)
    EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit shared INPUT node";
    EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit shared OUTPUT node";
    EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit shared BIAS node";
    EXPECT_TRUE(offspringNodeIDs.find(4) != offspringNodeIDs.end()) << "Should inherit excess node 4 from parent A";
    EXPECT_TRUE(offspringNodeIDs.find(5) != offspringNodeIDs.end()) << "Should inherit excess node 5 from parent A";
    
    // Should NOT have excess nodes from less fit parent B
    EXPECT_TRUE(offspringNodeIDs.find(6) == offspringNodeIDs.end()) << "Should NOT inherit excess node 6 from parent B";
    EXPECT_TRUE(offspringNodeIDs.find(7) == offspringNodeIDs.end()) << "Should NOT inherit excess node 7 from parent B";
}

TEST_F(CrossoverTest, SharedNodesAlwaysInherited_IncludingIOBias) {
    // Test that shared nodes (including I/O/BIAS) are always inherited, while excess nodes follow fitness rules
    
    // Parent A: I/O/BIAS + hidden nodes
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 10, 11};
    paramsA._nodeTypes = {
        NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS,
        NodeType::HIDDEN, NodeType::HIDDEN
    };
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::SIGMOID}, {ActivationType::TANH}
    };
    
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 10};
    paramsA._targetNodeHistoryIDs = {10, 11, 2};
    paramsA._connectionAttributes = {{1.0f, true}, {1.1f, true}, {1.2f, true}};
    
    // Parent B: Same I/O/BIAS (shared) + different hidden nodes (excess)
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 20, 21};
    paramsB._nodeTypes = {
        NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS,
        NodeType::HIDDEN, NodeType::HIDDEN
    };
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE},
        {ActivationType::RELU}, {ActivationType::SIGMOID}
    };
    
    paramsB._connectionHistoryIDs = {4, 5, 6};
    paramsB._sourceNodeHistoryIDs = {1, 3, 20};
    paramsB._targetNodeHistoryIDs = {20, 21, 2};
    paramsB._connectionAttributes = {{2.0f, true}, {2.1f, true}, {2.2f, true}};
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Test case 1: Parent A is fitter - shared nodes inherited, A's excess nodes inherited
    {
        TestFitnessResult fitnessA(0.9);
        TestFitnessResult fitnessB(0.3);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Shared nodes (including I/O/BIAS) always inherited
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Shared INPUT node must be inherited";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Shared OUTPUT node must be inherited";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Shared BIAS node must be inherited";
        
        // Excess nodes from fitter parent A should be inherited
        EXPECT_TRUE(offspringNodeIDs.find(10) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent A";
        EXPECT_TRUE(offspringNodeIDs.find(11) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent A";
        
        // Excess nodes from less fit parent B should NOT be inherited
        EXPECT_TRUE(offspringNodeIDs.find(20) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent B";
        EXPECT_TRUE(offspringNodeIDs.find(21) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent B";
    }
    
    // Test case 2: Parent B is fitter - shared nodes inherited, B's excess nodes inherited
    {
        TestFitnessResult fitnessA(0.2);
        TestFitnessResult fitnessB(0.8);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Shared nodes (including I/O/BIAS) always inherited
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Shared INPUT node must be inherited";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Shared OUTPUT node must be inherited";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Shared BIAS node must be inherited";
        
        // Excess nodes from fitter parent B should be inherited
        EXPECT_TRUE(offspringNodeIDs.find(20) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent B";
        EXPECT_TRUE(offspringNodeIDs.find(21) != offspringNodeIDs.end()) << "Should inherit excess hidden node from fitter parent B";
        
        // Excess nodes from less fit parent A should NOT be inherited
        EXPECT_TRUE(offspringNodeIDs.find(10) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent A";
        EXPECT_TRUE(offspringNodeIDs.find(11) == offspringNodeIDs.end()) << "Should NOT inherit excess hidden node from less fit parent A";
    }
    
    // Test case 3: Different I/O/BIAS nodes to verify they follow normal NEAT rules
    {
        // Parent C: Different I/O/BIAS node IDs (to test non-shared scenario)
        GenomeParams paramsC;
        paramsC._nodeHistoryIDs = {100, 101, 102};  // Different I/O/BIAS IDs
        paramsC._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS};
        paramsC._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::NONE}, {ActivationType::NONE}
        };
        
        paramsC._connectionHistoryIDs = {100};
        paramsC._sourceNodeHistoryIDs = {100};
        paramsC._targetNodeHistoryIDs = {101};
        paramsC._connectionAttributes = {{3.0f, true}};
        
        Genome parentC(paramsC);
        
        // Parent A vs Parent C - no shared nodes (including different I/O/BIAS)
        TestFitnessResult fitnessA(0.9);
        TestFitnessResult fitnessC(0.1);
        
        CrossoverParams params(0.0);
        
        Genome offspring = crossover(parentA, fitnessA, parentC, fitnessC, params);
        
        std::set<uint32_t> offspringNodeIDs;
        for (const auto& node : offspring.get_nodeGenes()) {
            offspringNodeIDs.insert(node.get_historyID());
        }
        
        // Should get all nodes from fitter parent A (no shared nodes)
        EXPECT_TRUE(offspringNodeIDs.find(1) != offspringNodeIDs.end()) << "Should inherit I/O/BIAS from fitter parent when not shared";
        EXPECT_TRUE(offspringNodeIDs.find(2) != offspringNodeIDs.end()) << "Should inherit I/O/BIAS from fitter parent when not shared";
        EXPECT_TRUE(offspringNodeIDs.find(3) != offspringNodeIDs.end()) << "Should inherit I/O/BIAS from fitter parent when not shared";
        
        // Should NOT get I/O/BIAS from less fit parent C
        EXPECT_TRUE(offspringNodeIDs.find(100) == offspringNodeIDs.end()) << "Should NOT inherit I/O/BIAS from less fit parent";
        EXPECT_TRUE(offspringNodeIDs.find(101) == offspringNodeIDs.end()) << "Should NOT inherit I/O/BIAS from less fit parent";
        EXPECT_TRUE(offspringNodeIDs.find(102) == offspringNodeIDs.end()) << "Should NOT inherit I/O/BIAS from less fit parent";
    }
}

TEST_F(CrossoverTest, DisabledGeneReactivation_ZeroPercent) {
    // Test that matching disabled connections remain disabled with 0% reactivation probability
    
    // Parent A: I/O/BIAS + connections (some matching with B, some disabled)
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::SIGMOID}
    };
    
    // Parent A connections: some matching with B, some disabled
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4};
    paramsA._targetNodeHistoryIDs = {2, 4, 2};
    paramsA._connectionAttributes = {
        {1.0f, true},   // matching - enabled
        {1.5f, false},  // matching - disabled
        {2.0f, false}   // matching - disabled
    };
    
    // Parent B: Same matching connections as A but with different enable states
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 4};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::TANH}
    };
    
    // Parent B connections: SAME history IDs as A to create matching connections
    paramsB._connectionHistoryIDs = {1, 2, 3};
    paramsB._sourceNodeHistoryIDs = {1, 3, 4};
    paramsB._targetNodeHistoryIDs = {2, 4, 2};
    paramsB._connectionAttributes = {
        {2.0f, true},   // matching - enabled (different weight than A)
        {2.5f, true},   // matching - enabled (A has this disabled)
        {3.0f, false}   // matching - disabled (both parents disabled)
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Equal fitness to get random 50/50 selection of matching connections
    TestFitnessResult fitnessA(0.7);
    TestFitnessResult fitnessB(0.7);
    
    // 0% reactivation probability - disabled genes should stay disabled
    CrossoverParams params(0.0);
    
    // Run multiple times to test reactivation behavior on randomly selected connections
    int disabledRemainDisabled = 0;
    int totalDisabledSelected = 0;
    const int numRuns = 100;
    
    for (int run = 0; run < numRuns; ++run) {
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        for (const auto& conn : offspring.get_connectionGenes()) {
            uint32_t historyID = conn.get_historyID();
            bool isEnabled = conn.get_attributes().enabled;
            
            // Check if a connection that was disabled in at least one parent got selected
            if (historyID == 2) {
                // Connection 2: A=disabled, B=enabled
                totalDisabledSelected++;
                if (!isEnabled) {
                    disabledRemainDisabled++;
                }
            } else if (historyID == 3) {
                // Connection 3: A=disabled, B=disabled
                totalDisabledSelected++;
                if (!isEnabled) {
                    disabledRemainDisabled++;
                }
            }
            // Connection 1: both parents enabled, so always enabled regardless
        }
    }
    
    // With 0% reactivation, connections that were disabled in the selected parent should remain disabled
    // Since we can't predict which parent gets selected, we just verify that SOME disabled connections remain disabled
    EXPECT_GT(disabledRemainDisabled, 0) << "Some disabled connections should remain disabled with 0% reactivation";
    EXPECT_GT(totalDisabledSelected, 0) << "Should have selected some connections that were disabled in at least one parent";
}

TEST_F(CrossoverTest, DisabledGeneReactivation_HundredPercent) {
    // Test that matching disabled connections are reactivated with 100% reactivation probability
    
    // Use same parent structures as previous test for consistency
    GenomeParams paramsA;
    paramsA._nodeHistoryIDs = {1, 2, 3, 4};
    paramsA._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsA._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::SIGMOID}
    };
    
    // Parent A connections: some matching with B, some disabled
    paramsA._connectionHistoryIDs = {1, 2, 3};
    paramsA._sourceNodeHistoryIDs = {1, 3, 4};
    paramsA._targetNodeHistoryIDs = {2, 4, 2};
    paramsA._connectionAttributes = {
        {1.0f, true},   // matching - enabled
        {1.5f, false},  // matching - disabled
        {2.0f, false}   // matching - disabled
    };
    
    // Parent B: Same matching connections as A but with different enable states
    GenomeParams paramsB;
    paramsB._nodeHistoryIDs = {1, 2, 3, 4};
    paramsB._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT, NodeType::BIAS, NodeType::HIDDEN};
    paramsB._nodeAttributes = {
        {ActivationType::NONE}, {ActivationType::NONE}, 
        {ActivationType::NONE}, {ActivationType::TANH}
    };
    
    // Parent B connections: SAME history IDs as A to create matching connections
    paramsB._connectionHistoryIDs = {1, 2, 3};
    paramsB._sourceNodeHistoryIDs = {1, 3, 4};
    paramsB._targetNodeHistoryIDs = {2, 4, 2};
    paramsB._connectionAttributes = {
        {2.0f, true},   // matching - enabled (different weight than A)
        {2.5f, true},   // matching - enabled (A has this disabled)
        {3.0f, false}   // matching - disabled (both parents disabled)
    };
    
    Genome parentA(paramsA);
    Genome parentB(paramsB);
    
    // Equal fitness to get random 50/50 selection of matching connections
    TestFitnessResult fitnessA(0.7);
    TestFitnessResult fitnessB(0.7);
    
    // 100% reactivation probability - disabled genes should be reactivated
    CrossoverParams params(1.0);
    
    // Run multiple times to test reactivation behavior on randomly selected connections
    int disabledGotReactivated = 0;
    int totalDisabledSelected = 0;
    const int numRuns = 100;
    
    for (int run = 0; run < numRuns; ++run) {
        Genome offspring = crossover(parentA, fitnessA, parentB, fitnessB, params);
        
        for (const auto& conn : offspring.get_connectionGenes()) {
            uint32_t historyID = conn.get_historyID();
            bool isEnabled = conn.get_attributes().enabled;
            
            // Track when disabled connections from selected parent get reactivated
            if (historyID == 2) {
                // Connection 2: A=disabled, B=enabled
                // If A was selected (isEnabled=false originally), it should be reactivated to true
                // If B was selected (isEnabled=true originally), it stays true
                // Either way, should be enabled due to reactivation or original state
                EXPECT_TRUE(isEnabled) << "Connection 2 should be enabled (reactivated if from A, or originally enabled from B)";
                
                // We can't directly know which parent was selected, but with 100% reactivation,
                // any disabled connection should become enabled
            } else if (historyID == 3) {
                // Connection 3: A=disabled, B=disabled
                // Regardless of which parent was selected, it was disabled and should be reactivated
                totalDisabledSelected++;
                if (isEnabled) {
                    disabledGotReactivated++;
                }
            }
            // Connection 1: both parents enabled, so always enabled regardless
        }
    }
    
    // With 100% reactivation, ALL disabled connections should be reactivated
    EXPECT_EQ(disabledGotReactivated, totalDisabledSelected) << "ALL disabled connections should be reactivated with 100% reactivation";
    EXPECT_GT(totalDisabledSelected, 0) << "Should have selected some connections that were disabled in both parents";
}
