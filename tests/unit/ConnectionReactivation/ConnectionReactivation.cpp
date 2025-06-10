#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/ConnectionReactivation.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"

using namespace Operator;

class ConnectionReactivationTest : public ::testing::Test {
protected:
    void SetUp() override {
        std::srand(42); // Fixed seed for reproducibility where needed
    }

    // Helper to create genome with some disabled connections
    Genome createGenomeWithDisabledConnections() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3, 4};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, 
            {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1, 2, 3, 4};
        params._sourceNodeHistoryIDs = {1, 1, 2, 3};
        params._targetNodeHistoryIDs = {2, 3, 3, 4};
        params._connectionAttributes = {
            {1.0f, true},   // enabled
            {1.5f, false},  // disabled
            {2.0f, false},  // disabled  
            {2.5f, true}    // enabled
        };
        return Genome(params);
    }

    // Helper to create genome with all connections enabled
    Genome createGenomeAllEnabled() {
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
            {1.0f, true},   // enabled
            {1.5f, true}    // enabled
        };
        return Genome(params);
    }

    // Helper to create genome with single disabled connection
    Genome createGenomeSingleDisabled() {
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
            {1.0f, true},   // enabled
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

    // Helper to verify genome structural integrity
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
        
        // Same connection structure (except enabled flags)
        for (size_t i = 0; i < original.get_connectionGenes().size(); ++i) {
            const auto& origConn = original.get_connectionGenes()[i];
            const auto& mutConn = mutated.get_connectionGenes()[i];
            EXPECT_EQ(origConn.get_historyID(), mutConn.get_historyID());
            EXPECT_EQ(origConn.get_sourceNodeGene().get_historyID(), 
                     mutConn.get_sourceNodeGene().get_historyID());
            EXPECT_EQ(origConn.get_targetNodeGene().get_historyID(), 
                     mutConn.get_targetNodeGene().get_historyID());
            EXPECT_FLOAT_EQ(origConn.get_attributes().weight, mutConn.get_attributes().weight);
        }
    }
};

// ============================================================================
// Basic Reactivation Mechanics Tests
// ============================================================================

TEST_F(ConnectionReactivationTest, ConnectionAlwaysReactivated) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    size_t originalEnabled = countEnabledConnections(original);
    size_t originalDisabled = countDisabledConnections(original);
    
    Genome mutated = connectionReactivation(original, params);
    
    // Should have exactly one more enabled connection
    EXPECT_EQ(countEnabledConnections(mutated), originalEnabled + 1);
    EXPECT_EQ(countDisabledConnections(mutated), originalDisabled - 1);
}

TEST_F(ConnectionReactivationTest, ReactivationCountIncrement_Deterministic) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    size_t originalEnabled = countEnabledConnections(original);
    
    // Run multiple times - should always reactivate exactly one connection
    for (int i = 0; i < 10; ++i) {
        Genome mutated = connectionReactivation(original, params);
        EXPECT_EQ(countEnabledConnections(mutated), originalEnabled + 1);
    }
}

TEST_F(ConnectionReactivationTest, NoStructuralChanges) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    Genome mutated = connectionReactivation(original, params);
    
    // Same number of connections, same innovation numbers
    EXPECT_EQ(original.get_connectionGenes().size(), mutated.get_connectionGenes().size());
    
    for (size_t i = 0; i < original.get_connectionGenes().size(); ++i) {
        EXPECT_EQ(original.get_connectionGenes()[i].get_historyID(), 
                 mutated.get_connectionGenes()[i].get_historyID());
    }
}

// ============================================================================
// Selection Strategy Tests
// ============================================================================

TEST_F(ConnectionReactivationTest, RandomStrategy_CanProduceDifferent) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    std::set<uint32_t> reactivatedInnovations;
    
    // Run enough iterations to likely see different selections
    for (int i = 0; i < 50; ++i) {
        Genome mutated = connectionReactivation(original, params);
        
        // Find which connection was reactivated
        for (size_t j = 0; j < original.get_connectionGenes().size(); ++j) {
            const auto& origConn = original.get_connectionGenes()[j];
            const auto& mutConn = mutated.get_connectionGenes()[j];
            
            if (!origConn.get_attributes().enabled && mutConn.get_attributes().enabled) {
                reactivatedInnovations.insert(origConn.get_historyID());
                break;
            }
        }
    }
    
    // Should observe more than one unique reactivation (accounting for randomness outliers)
    EXPECT_GT(reactivatedInnovations.size(), 1u);
}

TEST_F(ConnectionReactivationTest, OldestFirstStrategy_Deterministic) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::OLDEST_FIRST);
    
    // Find the lowest innovation number among disabled connections
    uint32_t expectedInnovation = UINT32_MAX;
    for (const auto& conn : original.get_connectionGenes()) {
        if (!conn.get_attributes().enabled) {
            expectedInnovation = std::min(expectedInnovation, conn.get_historyID());
        }
    }
    
    // Run multiple times - should always reactivate the same (oldest) connection
    for (int i = 0; i < 10; ++i) {
        Genome mutated = connectionReactivation(original, params);
        
        // Find which connection was reactivated
        bool found = false;
        for (size_t j = 0; j < original.get_connectionGenes().size(); ++j) {
            const auto& origConn = original.get_connectionGenes()[j];
            const auto& mutConn = mutated.get_connectionGenes()[j];
            
            if (!origConn.get_attributes().enabled && mutConn.get_attributes().enabled) {
                EXPECT_EQ(origConn.get_historyID(), expectedInnovation);
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
}

TEST_F(ConnectionReactivationTest, NewestFirstStrategy_Deterministic) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::NEWEST_FIRST);
    
    // Find the highest innovation number among disabled connections
    uint32_t expectedInnovation = 0;
    for (const auto& conn : original.get_connectionGenes()) {
        if (!conn.get_attributes().enabled) {
            expectedInnovation = std::max(expectedInnovation, conn.get_historyID());
        }
    }
    
    // Run multiple times - should always reactivate the same (newest) connection
    for (int i = 0; i < 10; ++i) {
        Genome mutated = connectionReactivation(original, params);
        
        // Find which connection was reactivated
        bool found = false;
        for (size_t j = 0; j < original.get_connectionGenes().size(); ++j) {
            const auto& origConn = original.get_connectionGenes()[j];
            const auto& mutConn = mutated.get_connectionGenes()[j];
            
            if (!origConn.get_attributes().enabled && mutConn.get_attributes().enabled) {
                EXPECT_EQ(origConn.get_historyID(), expectedInnovation);
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
}

TEST_F(ConnectionReactivationTest, SingleDisabledConnection_Deterministic) {
    Genome original = createGenomeSingleDisabled();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    // Find the disabled connection
    uint32_t disabledInnovation = 0;
    for (const auto& conn : original.get_connectionGenes()) {
        if (!conn.get_attributes().enabled) {
            disabledInnovation = conn.get_historyID();
            break;
        }
    }
    
    // Run multiple times - should always reactivate the same connection
    for (int i = 0; i < 10; ++i) {
        Genome mutated = connectionReactivation(original, params);
        
        // Verify the specific connection was reactivated
        bool found = false;
        for (size_t j = 0; j < original.get_connectionGenes().size(); ++j) {
            const auto& origConn = original.get_connectionGenes()[j];
            const auto& mutConn = mutated.get_connectionGenes()[j];
            
            if (!origConn.get_attributes().enabled && mutConn.get_attributes().enabled) {
                EXPECT_EQ(origConn.get_historyID(), disabledInnovation);
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
}

// ============================================================================
// Genome Integrity Tests
// ============================================================================

TEST_F(ConnectionReactivationTest, OriginalGenomeUntouched) {
    Genome original = createGenomeWithDisabledConnections();
    
    // Store original enabled counts
    size_t origEnabled = countEnabledConnections(original);
    size_t origDisabled = countDisabledConnections(original);
    
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    Genome mutated = connectionReactivation(original, params);
    
    // Verify original unchanged
    EXPECT_EQ(countEnabledConnections(original), origEnabled);
    EXPECT_EQ(countDisabledConnections(original), origDisabled);
}

TEST_F(ConnectionReactivationTest, AllExistingDataPreserved) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    Genome mutated = connectionReactivation(original, params);
    
    verifyGenomeIntegrity(original, mutated);
}

TEST_F(ConnectionReactivationTest, OnlyEnabledFlagChanges) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    Genome mutated = connectionReactivation(original, params);
    
    // Verify only one connection changed enabled state, everything else identical
    int changedCount = 0;
    for (size_t i = 0; i < original.get_connectionGenes().size(); ++i) {
        const auto& origConn = original.get_connectionGenes()[i];
        const auto& mutConn = mutated.get_connectionGenes()[i];
        
        // Same innovation, nodes, weight
        EXPECT_EQ(origConn.get_historyID(), mutConn.get_historyID());
        EXPECT_FLOAT_EQ(origConn.get_attributes().weight, mutConn.get_attributes().weight);
        
        // Count enabled state changes
        if (origConn.get_attributes().enabled != mutConn.get_attributes().enabled) {
            changedCount++;
            // Should change from disabled to enabled
            EXPECT_FALSE(origConn.get_attributes().enabled);
            EXPECT_TRUE(mutConn.get_attributes().enabled);
        }
    }
    
    // Exactly one connection should have changed
    EXPECT_EQ(changedCount, 1);
}

TEST_F(ConnectionReactivationTest, ValidConstruction) {
    Genome original = createGenomeWithDisabledConnections();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    Genome mutated = connectionReactivation(original, params);
    
    // Should be able to construct phenotype without errors
    EXPECT_NO_THROW({
        Operator::phenotypeConstruct(mutated);
        const auto& phenotype = mutated.get_phenotype();
        EXPECT_FALSE(phenotype._nodeGeneAttributes.empty());
    });
}

// ============================================================================
// Failure Cases (User Error Assertions) Tests
// ============================================================================

TEST_F(ConnectionReactivationTest, NoDisabledConnections_AllEnabled) {
    Genome original = createGenomeAllEnabled();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    // Should assert when no disabled connections available
    EXPECT_DEATH({
        connectionReactivation(original, params);
    }, "");
}

TEST_F(ConnectionReactivationTest, NoConnections_EmptyGenome) {
    Genome original = createEmptyGenome();
    ConnectionReactivationParams params(ConnectionReactivationParams::SelectionStrategy::RANDOM);
    
    // Should assert when no connections exist
    EXPECT_DEATH({
        connectionReactivation(original, params);
    }, "");
}

// ============================================================================
// Strategy Parameter Tests
// ============================================================================

TEST_F(ConnectionReactivationTest, ValidStrategies) {
    EXPECT_NO_THROW({
        ConnectionReactivationParams params1(ConnectionReactivationParams::SelectionStrategy::RANDOM);
        ConnectionReactivationParams params2(ConnectionReactivationParams::SelectionStrategy::OLDEST_FIRST);
        ConnectionReactivationParams params3(ConnectionReactivationParams::SelectionStrategy::NEWEST_FIRST);
    });
}