#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include <unordered_set>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/ConnectionMutation.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Operator;

class ConnectionMutationTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_unique<HistoryTracker>();
        std::srand(42); // Fixed seed for reproducibility where needed
    }

    std::unique_ptr<HistoryTracker> historyTracker;

    // Helper to create basic genome with disconnected nodes
    Genome createDisconnectedGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3, 4};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, 
            {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        // No connections
        return Genome(params);
    }

    // Helper to create genome with some connections
    Genome createPartiallyConnectedGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3, 4};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, 
            {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        params._connectionHistoryIDs = {1};
        params._sourceNodeHistoryIDs = {1};
        params._targetNodeHistoryIDs = {2};
        params._connectionAttributes = {{1.0f, true}};
        return Genome(params);
    }

    // Helper to create fully connected genome (all possible connections)
    Genome createFullyConnectedGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1, 2, 3};
        params._nodeTypes = {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT};
        params._nodeAttributes = {
            {ActivationType::NONE}, {ActivationType::SIGMOID}, {ActivationType::NONE}
        };
        
        // Add all possible connections (excluding self and invalid combinations)
        // 1->2, 1->3, 2->3
        params._connectionHistoryIDs = {1, 2, 3};
        params._sourceNodeHistoryIDs = {1, 1, 2};
        params._targetNodeHistoryIDs = {2, 3, 3};
        params._connectionAttributes = {
            {1.0f, true}, {1.5f, true}, {2.0f, true}
        };
        return Genome(params);
    }

    // Helper to create single node genome
    Genome createSingleNodeGenome() {
        GenomeParams params;
        params._nodeHistoryIDs = {1};
        params._nodeTypes = {NodeType::HIDDEN};
        params._nodeAttributes = {{ActivationType::SIGMOID}};
        return Genome(params);
    }

    // Helper to create empty genome
    Genome createEmptyGenome() {
        GenomeParams params;
        return Genome(params);
    }

    // Helper to count unique connection pairs in a genome
    std::set<std::pair<uint32_t, uint32_t>> getConnectionPairs(const Genome& genome) {
        std::set<std::pair<uint32_t, uint32_t>> pairs;
        for (const auto& conn : genome.get_connectionGenes()) {
            pairs.insert({conn.get_sourceNodeGene().get_historyID(), 
                         conn.get_targetNodeGene().get_historyID()});
        }
        return pairs;
    }

    // Helper to verify genome structural integrity
    void verifyGenomeIntegrity(const Genome& original, const Genome& mutated) {
        // Same number of nodes
        EXPECT_EQ(original.get_nodeGenes().size(), mutated.get_nodeGenes().size());
        
        // Same node data
        for (size_t i = 0; i < original.get_nodeGenes().size(); ++i) {
            const auto& origNode = original.get_nodeGenes()[i];
            const auto& mutNode = mutated.get_nodeGenes()[i];
            EXPECT_EQ(origNode.get_historyID(), mutNode.get_historyID());
            EXPECT_EQ(origNode.get_type(), mutNode.get_type());
            EXPECT_EQ(origNode.get_attributes().activationType, mutNode.get_attributes().activationType);
        }
        
        // Original connections preserved
        const auto& origConns = original.get_connectionGenes();
        const auto& mutConns = mutated.get_connectionGenes();
        
        for (size_t i = 0; i < origConns.size(); ++i) {
            const auto& origConn = origConns[i];
            // Find matching connection in mutated genome
            bool found = false;
            for (const auto& mutConn : mutConns) {
                if (origConn.get_historyID() == mutConn.get_historyID()) {
                    EXPECT_EQ(origConn.get_sourceNodeGene().get_historyID(), 
                             mutConn.get_sourceNodeGene().get_historyID());
                    EXPECT_EQ(origConn.get_targetNodeGene().get_historyID(), 
                             mutConn.get_targetNodeGene().get_historyID());
                    EXPECT_FLOAT_EQ(origConn.get_attributes().weight, mutConn.get_attributes().weight);
                    EXPECT_EQ(origConn.get_attributes().enabled, mutConn.get_attributes().enabled);
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "Original connection " << origConn.get_historyID() << " not found in mutated genome";
        }
    }
};

// ============================================================================
// Basic Creation Mechanics Tests
// ============================================================================

TEST_F(ConnectionMutationTest, ConnectionAlwaysCreated_DisconnectedGenome) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // Should have exactly one more connection
    EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 1);
}

TEST_F(ConnectionMutationTest, ConnectionAlwaysCreated_PartiallyConnected) {
    Genome original = createPartiallyConnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::RECURRENT);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // Should have exactly one more connection
    EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 1);
}

TEST_F(ConnectionMutationTest, ConnectionCountIncrement_Deterministic) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(0.5, 1.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Run multiple times - should always add exactly one connection
    for (int i = 0; i < 10; ++i) {
        auto localHistoryTracker = std::make_unique<HistoryTracker>();
        Genome mutated = connectionMutation(original, std::move(localHistoryTracker), params);
        EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 1);
    }
}

TEST_F(ConnectionMutationTest, FreshInnovationNumber) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // New connection should have a valid innovation number
    EXPECT_EQ(mutated.get_connectionGenes().size(), 1u);
    const auto& newConn = mutated.get_connectionGenes()[0];
    EXPECT_GT(newConn.get_historyID(), 0u);
}

// ============================================================================
// Connection Properties Tests
// ============================================================================

TEST_F(ConnectionMutationTest, WeightWithinRange) {
    Genome original = createDisconnectedGenome();
    double weightRange = 1.5;
    ConnectionMutationParams params(1.0, weightRange, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Test multiple times to check range compliance
    for (int i = 0; i < 50; ++i) {
        historyTracker = std::make_unique<HistoryTracker>(); // Reset for fresh innovation numbers
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
        
        const auto& newConn = mutated.get_connectionGenes()[0];
        float weight = newConn.get_attributes().weight;
        EXPECT_GE(weight, -weightRange);
        EXPECT_LE(weight, weightRange);
    }
}

TEST_F(ConnectionMutationTest, AlwaysEnabled) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    for (int i = 0; i < 10; ++i) {
        historyTracker = std::make_unique<HistoryTracker>();
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
        
        const auto& newConn = mutated.get_connectionGenes()[0];
        EXPECT_TRUE(newConn.get_attributes().enabled);
    }
}

TEST_F(ConnectionMutationTest, NoSelfConnections) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::RECURRENT);
    
    // Run many times to verify no self-connections are created
    for (int i = 0; i < 100; ++i) {
        historyTracker = std::make_unique<HistoryTracker>();
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
        
        const auto& newConn = mutated.get_connectionGenes()[0];
        EXPECT_NE(newConn.get_sourceNodeGene().get_historyID(), 
                 newConn.get_targetNodeGene().get_historyID());
    }
}

TEST_F(ConnectionMutationTest, NoDuplicateConnections) {
    Genome original = createPartiallyConnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    auto originalPairs = getConnectionPairs(original);
    
    for (int i = 0; i < 20; ++i) {
        historyTracker = std::make_unique<HistoryTracker>();
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
        
        auto mutatedPairs = getConnectionPairs(mutated);
        
        // Should have one new pair that wasn't in original
        EXPECT_EQ(mutatedPairs.size(), originalPairs.size() + 1);
        
        // All original pairs should still exist
        for (const auto& origPair : originalPairs) {
            EXPECT_TRUE(mutatedPairs.find(origPair) != mutatedPairs.end());
        }
    }
}

// ============================================================================
// Randomness Testing Tests
// ============================================================================

TEST_F(ConnectionMutationTest, MultiplePossibleOutcomes_CanProduceDifferent) {
    Genome original = createDisconnectedGenome(); // 4 nodes, many possible connections
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::RECURRENT);
    
    std::set<std::pair<uint32_t, uint32_t>> observedPairs;
    
    // Run enough iterations to likely see different pairs
    for (int i = 0; i < 100; ++i) {
        historyTracker = std::make_unique<HistoryTracker>();
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
        
        const auto& newConn = mutated.get_connectionGenes()[0];
        observedPairs.insert({newConn.get_sourceNodeGene().get_historyID(),
                             newConn.get_targetNodeGene().get_historyID()});
    }
    
    // Should observe more than one unique pair (accounting for randomness outliers)
    EXPECT_GT(observedPairs.size(), 1u);
}

TEST_F(ConnectionMutationTest, DeterministicWhenSingleChoice) {
    // Create genome where only one connection was possible, but it already exists
    GenomeParams params;
    params._nodeHistoryIDs = {1, 2};
    params._nodeTypes = {NodeType::INPUT, NodeType::OUTPUT};
    params._nodeAttributes = {{ActivationType::NONE}, {ActivationType::NONE}};
    params._connectionHistoryIDs = {1};
    params._sourceNodeHistoryIDs = {1};
    params._targetNodeHistoryIDs = {2};
    params._connectionAttributes = {{1.0f, true}};
    Genome original(params);
    
    ConnectionMutationParams mutParams(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Should assert because no more connections are possible
    EXPECT_DEATH({
        connectionMutation(original, std::move(historyTracker), mutParams);
    }, "");
}

// ============================================================================
// Genome Integrity Tests
// ============================================================================

TEST_F(ConnectionMutationTest, OriginalGenomeUntouched) {
    Genome original = createPartiallyConnectedGenome();
    
    // Store original data
    size_t origNodeCount = original.get_nodeGenes().size();
    size_t origConnCount = original.get_connectionGenes().size();
    auto origFirstConnWeight = original.get_connectionGenes()[0].get_attributes().weight;
    
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // Verify original unchanged
    EXPECT_EQ(original.get_nodeGenes().size(), origNodeCount);
    EXPECT_EQ(original.get_connectionGenes().size(), origConnCount);
    EXPECT_FLOAT_EQ(original.get_connectionGenes()[0].get_attributes().weight, origFirstConnWeight);
}

TEST_F(ConnectionMutationTest, AllExistingDataPreserved) {
    Genome original = createPartiallyConnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    verifyGenomeIntegrity(original, mutated);
}

TEST_F(ConnectionMutationTest, OnlyAdditionOccurs) {
    Genome original = createPartiallyConnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // Exactly one more connection, everything else identical
    EXPECT_EQ(mutated.get_nodeGenes().size(), original.get_nodeGenes().size());
    EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 1);
}

TEST_F(ConnectionMutationTest, ValidConstruction) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    
    // Should be able to construct phenotype without errors
    EXPECT_NO_THROW({
        Operator::phenotypeConstruct(mutated);
        const auto& phenotype = mutated.get_phenotype();
        EXPECT_FALSE(phenotype._nodeGeneAttributes.empty());
    });
}

// ============================================================================
// Parameter Validation Tests
// ============================================================================

TEST_F(ConnectionMutationTest, InvalidRateRejection_Zero) {
    EXPECT_DEATH({
        ConnectionMutationParams params(0.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    }, "");
}

TEST_F(ConnectionMutationTest, InvalidRateRejection_Negative) {
    EXPECT_DEATH({
        ConnectionMutationParams params(-0.1, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    }, "");
}

TEST_F(ConnectionMutationTest, InvalidRateRejection_TooHigh) {
    EXPECT_DEATH({
        ConnectionMutationParams params(1.1, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    }, "");
}

TEST_F(ConnectionMutationTest, InvalidWeightRangeRejection_Zero) {
    EXPECT_DEATH({
        ConnectionMutationParams params(0.5, 0.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    }, "");
}

TEST_F(ConnectionMutationTest, InvalidWeightRangeRejection_Negative) {
    EXPECT_DEATH({
        ConnectionMutationParams params(0.5, -1.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    }, "");
}

TEST_F(ConnectionMutationTest, ValidBoundaryValues) {
    EXPECT_NO_THROW({
        ConnectionMutationParams params1(0.001, 0.001, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
        ConnectionMutationParams params2(1.0, 100.0, ConnectionMutationParams::NetworkTopology::RECURRENT);
    });
}

// ============================================================================
// Failure Cases (User Error Assertions) Tests
// ============================================================================

TEST_F(ConnectionMutationTest, InsufficientNodes_Empty) {
    Genome original = createEmptyGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Should assert with 0 nodes
    EXPECT_DEATH({
        connectionMutation(original, std::move(historyTracker), params);
    }, "");
}

TEST_F(ConnectionMutationTest, InsufficientNodes_Single) {
    Genome original = createSingleNodeGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Should assert with single node
    EXPECT_DEATH({
        connectionMutation(original, std::move(historyTracker), params);
    }, "");
}

TEST_F(ConnectionMutationTest, NoAvailablePairs_FullyConnected) {
    Genome original = createFullyConnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Should assert when no more connections possible
    EXPECT_DEATH({
        connectionMutation(original, std::move(historyTracker), params);
    }, "");
}

TEST_F(ConnectionMutationTest, ValidHistoryTracker) {
    Genome original = createDisconnectedGenome();
    ConnectionMutationParams params(1.0, 2.0, ConnectionMutationParams::NetworkTopology::FEED_FORWARD);
    
    // Test with valid history tracker (should work)
    EXPECT_NO_THROW({
        Genome mutated = connectionMutation(original, std::move(historyTracker), params);
    });
}

// ============================================================================
// Large Genome Performance Test
// ============================================================================

TEST_F(ConnectionMutationTest, LargeGenome_Performance) {
    // Create large genome with many nodes, some connections
    GenomeParams params;
    const size_t numNodes = 50;
    
    for (size_t i = 1; i <= numNodes; ++i) {
        params._nodeHistoryIDs.push_back(static_cast<uint32_t>(i));
        if (i <= 10) {
            params._nodeTypes.push_back(NodeType::INPUT);
        } else if (i > 40) {
            params._nodeTypes.push_back(NodeType::OUTPUT);
        } else {
            params._nodeTypes.push_back(NodeType::HIDDEN);
        }
        params._nodeAttributes.push_back({ActivationType::SIGMOID});
    }
    
    // Add some existing connections
    for (size_t i = 0; i < 20; ++i) {
        params._connectionHistoryIDs.push_back(static_cast<uint32_t>(i + 1));
        params._sourceNodeHistoryIDs.push_back(static_cast<uint32_t>(i + 1));
        params._targetNodeHistoryIDs.push_back(static_cast<uint32_t>(i + 11));
        params._connectionAttributes.push_back({1.0f, true});
    }
    
    Genome original(params);
    ConnectionMutationParams mutParams(1.0, 2.0, ConnectionMutationParams::NetworkTopology::RECURRENT);
    
    // Should handle large genome efficiently
    EXPECT_NO_THROW({
        Genome mutated = connectionMutation(original, std::move(historyTracker), mutParams);
        EXPECT_EQ(mutated.get_connectionGenes().size(), original.get_connectionGenes().size() + 1);
    });
}