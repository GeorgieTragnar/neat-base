#include <gtest/gtest.h>
#include "core/Genome.hpp"
#include "core/Gene.hpp"

namespace {

class GenomeTest : public ::testing::Test {
protected:
    neat::core::Genome::Config config;
    
    void SetUp() override {
        config = neat::core::Genome::Config(2, 1); // 2 inputs, 1 output
    }
};

TEST_F(GenomeTest, InitialTopology) {
    // Arrange
    neat::core::Genome genome(config);
    
    // Act
    genome.initMinimalTopology(2, 1);
    
    // Assert
    const auto& nodes = genome.getNodes();
    const auto& genes = genome.getGenes();
    
    // Check nodes
    EXPECT_EQ(4, nodes.size()); // 2 inputs + 1 output + 1 bias
    
    int inputCount = 0, outputCount = 0, biasCount = 0;
    for (const auto& [id, type] : nodes) {
        if (type == neat::core::ENodeType::INPUT) inputCount++;
        if (type == neat::core::ENodeType::OUTPUT) outputCount++;
        if (type == neat::core::ENodeType::BIAS) biasCount++;
    }
    
    EXPECT_EQ(2, inputCount);
    EXPECT_EQ(1, outputCount);
    EXPECT_EQ(1, biasCount);
    
    // Check genes - should have connections from all inputs and bias to output
    EXPECT_EQ(3, genes.size()); // 2 inputs + 1 bias to 1 output
    
    // All genes should be enabled
    for (const auto& gene : genes) {
        EXPECT_TRUE(gene.enabled);
    }
}

TEST_F(GenomeTest, Validation) {
    // Arrange
    neat::core::Genome genome(config);
    genome.initMinimalTopology(2, 1);
    
    // Act & Assert
    EXPECT_NO_THROW(genome.validate());
}

// Add tests for mutation operations
TEST_F(GenomeTest, AddNodeMutation) {
    // Arrange
    neat::core::Genome genome(config);
    genome.initMinimalTopology(2, 1);
    size_t initialNodeCount = genome.getNodes().size();
    size_t initialGeneCount = genome.getGenes().size();
    
    // Act
    bool mutationApplied = genome.addNodeMutation();
    
    // Assert
    EXPECT_TRUE(mutationApplied);
    EXPECT_EQ(initialNodeCount + 1, genome.getNodes().size());
    // Fixed: addNodeMutation adds two new genes (keeping the original disabled)
    EXPECT_EQ(initialGeneCount + 2, genome.getGenes().size()); 
    
    // One gene should be disabled now
    int disabledCount = 0;
    for (const auto& gene : genome.getGenes()) {
        if (!gene.enabled) disabledCount++;
    }
    EXPECT_EQ(1, disabledCount);
    
    // Should still be valid
    EXPECT_NO_THROW(genome.validate());
}

} // namespace