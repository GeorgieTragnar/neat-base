#include <gtest/gtest.h>
#include <memory>
#include <cstdlib>
#include <ctime>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/operator/Init.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Operator;

class InitTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_shared<HistoryTracker>();
        std::srand(static_cast<unsigned>(std::time(nullptr)));
    }

    std::shared_ptr<HistoryTracker> historyTracker;
    
    std::vector<NodeGeneAttributes> createInputAttributes(size_t count) {
        std::vector<NodeGeneAttributes> attrs;
        for (size_t i = 0; i < count; i++) {
            attrs.push_back(NodeGeneAttributes{ActivationType::NONE});
        }
        return attrs;
    }
    
    std::vector<NodeGeneAttributes> createOutputAttributes(size_t count) {
        std::vector<NodeGeneAttributes> attrs;
        for (size_t i = 0; i < count; i++) {
            attrs.push_back(NodeGeneAttributes{ActivationType::SIGMOID});
        }
        return attrs;
    }
};

TEST_F(InitTest, CONNECT_TO_OUTPUTS_BasicFunctionality) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    EXPECT_EQ(genome.get_nodeGenes().size(), 5u); // 2 inputs + 1 bias + 2 outputs
    EXPECT_EQ(genome.get_connectionGenes().size(), 6u); // 2x2 input-output + 2 bias-output
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_AllInputsToAllOutputs) {
    auto inputAttrs = createInputAttributes(3);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    size_t expectedConnections = (3 * 2) + 2; // inputs×outputs + bias×outputs
    EXPECT_EQ(genome.get_connectionGenes().size(), expectedConnections);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_CorrectNodeTypes) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& nodes = genome.get_nodeGenes();
    
    // First 2 should be inputs
    EXPECT_EQ(nodes[0].get_type(), NodeType::INPUT);
    EXPECT_EQ(nodes[1].get_type(), NodeType::INPUT);
    
    // Third should be bias
    EXPECT_EQ(nodes[2].get_type(), NodeType::BIAS);
    
    // Last 2 should be outputs
    EXPECT_EQ(nodes[3].get_type(), NodeType::OUTPUT);
    EXPECT_EQ(nodes[4].get_type(), NodeType::OUTPUT);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_CorrectNodeAttributes) {
    auto inputAttrs = createInputAttributes(1);
    inputAttrs[0].activationType = ActivationType::RELU;
    
    auto outputAttrs = createOutputAttributes(1);
    outputAttrs[0].activationType = ActivationType::TANH;
    
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& nodes = genome.get_nodeGenes();
    
    EXPECT_EQ(nodes[0].get_attributes().activationType, ActivationType::RELU);
    EXPECT_EQ(nodes[1].get_attributes().activationType, ActivationType::NONE); // Bias
    EXPECT_EQ(nodes[2].get_attributes().activationType, ActivationType::TANH);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_AllConnectionsEnabled) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    for (const auto& conn : connections) {
        EXPECT_TRUE(conn.get_attributes().enabled);
    }
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_WeightsInRange) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    for (size_t i = 0; i < 4; i++) { // Only check input-output connections, not bias
        float weight = connections[i].get_attributes().weight;
        EXPECT_GE(weight, -2.0f);
        EXPECT_LE(weight, 2.0f);
    }
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_BiasDefaultAttributes) {
    auto inputAttrs = createInputAttributes(1);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs; // Empty - use defaults
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    // Last 2 connections should be bias connections with default weight 1.0
    EXPECT_EQ(connections[connections.size()-2].get_attributes().weight, 1.0f);
    EXPECT_EQ(connections[connections.size()-1].get_attributes().weight, 1.0f);
    EXPECT_TRUE(connections[connections.size()-2].get_attributes().enabled);
    EXPECT_TRUE(connections[connections.size()-1].get_attributes().enabled);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_BiasCustomAttributes) {
    auto inputAttrs = createInputAttributes(1);
    auto outputAttrs = createOutputAttributes(2);
    
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    biasAttrs[0] = ConnectionGeneAttributes{2.5f, false}; // Custom weight and disabled
    biasAttrs[1] = ConnectionGeneAttributes{-1.5f, true}; // Custom weight and enabled
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    // Last 2 connections should be bias connections with custom attributes
    EXPECT_EQ(connections[connections.size()-2].get_attributes().weight, 2.5f);
    EXPECT_FALSE(connections[connections.size()-2].get_attributes().enabled);
    EXPECT_EQ(connections[connections.size()-1].get_attributes().weight, -1.5f);
    EXPECT_TRUE(connections[connections.size()-1].get_attributes().enabled);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_BiasPartialAttributes) {
    auto inputAttrs = createInputAttributes(1);
    auto outputAttrs = createOutputAttributes(3);
    
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    biasAttrs[1] = ConnectionGeneAttributes{3.0f, false}; // Only middle output has custom bias
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    // Last 3 connections should be bias connections
    size_t lastIdx = connections.size() - 1;
    EXPECT_EQ(connections[lastIdx-2].get_attributes().weight, 1.0f); // Default
    EXPECT_TRUE(connections[lastIdx-2].get_attributes().enabled);
    
    EXPECT_EQ(connections[lastIdx-1].get_attributes().weight, 3.0f); // Custom
    EXPECT_FALSE(connections[lastIdx-1].get_attributes().enabled);
    
    EXPECT_EQ(connections[lastIdx].get_attributes().weight, 1.0f); // Default
    EXPECT_TRUE(connections[lastIdx].get_attributes().enabled);
}

TEST_F(InitTest, CONNECT_TO_OUTPUTS_UniqueInnovationNumbers) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    const auto& connections = genome.get_connectionGenes();
    
    std::vector<uint32_t> innovationNumbers;
    for (const auto& conn : connections) {
        innovationNumbers.push_back(conn.get_historyID());
    }
    
    std::sort(innovationNumbers.begin(), innovationNumbers.end());
    auto last = std::unique(innovationNumbers.begin(), innovationNumbers.end());
    EXPECT_EQ(last, innovationNumbers.end()); // No duplicates
}

TEST_F(InitTest, NONE_Strategy_NoConnections) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::NONE);
    
    Genome genome = init(historyTracker, params);
    
    EXPECT_EQ(genome.get_nodeGenes().size(), 5u); // Nodes still created
    EXPECT_EQ(genome.get_connectionGenes().size(), 0u); // No connections
}

TEST_F(InitTest, EdgeCase_SingleInputSingleOutput) {
    auto inputAttrs = createInputAttributes(1);
    auto outputAttrs = createOutputAttributes(1);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    EXPECT_EQ(genome.get_nodeGenes().size(), 3u); // 1 input + 1 bias + 1 output
    EXPECT_EQ(genome.get_connectionGenes().size(), 2u); // 1 input-output + 1 bias-output
}

TEST_F(InitTest, EdgeCase_EmptyInputs) {
    auto inputAttrs = createInputAttributes(0);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    EXPECT_EQ(genome.get_nodeGenes().size(), 3u); // 0 inputs + 1 bias + 2 outputs
    EXPECT_EQ(genome.get_connectionGenes().size(), 2u); // Only bias-output connections
}

TEST_F(InitTest, EdgeCase_EmptyOutputs) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(0);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    EXPECT_EQ(genome.get_nodeGenes().size(), 3u); // 2 inputs + 1 bias + 0 outputs
    EXPECT_EQ(genome.get_connectionGenes().size(), 0u); // No connections possible
}

TEST_F(InitTest, PhenotypeConstructionReadiness) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    Genome genome = init(historyTracker, params);
    
    // Should be able to construct phenotype without throwing
    EXPECT_NO_THROW({
        Operator::phenotypeConstruct(genome);
        const auto& phenotype = genome.get_phenotype();
    });
}

TEST_F(InitTest, BiasAttributes_ExtraIndicesIgnored) {
    auto inputAttrs = createInputAttributes(1);
    auto outputAttrs = createOutputAttributes(2);
    
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    biasAttrs[0] = ConnectionGeneAttributes{2.0f, true};
    biasAttrs[5] = ConnectionGeneAttributes{99.0f, false}; // Non-existent output index
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    // Should not throw despite extra bias attribute
    EXPECT_NO_THROW({
        Genome genome = init(historyTracker, params);
    });
}

TEST_F(InitTest, ConsistentHistoryIDs) {
    auto inputAttrs = createInputAttributes(2);
    auto outputAttrs = createOutputAttributes(2);
    std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
    
    InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                               InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
    
    // Create two genomes with same parameters using same history tracker
    Genome genome1 = init(historyTracker, params);
    
    // Reset history tracker to test consistency
    historyTracker = std::make_shared<HistoryTracker>();
    Genome genome2 = init(historyTracker, params);
    
    // Should have identical structure
    EXPECT_EQ(genome1.get_nodeGenes().size(), genome2.get_nodeGenes().size());
    EXPECT_EQ(genome1.get_connectionGenes().size(), genome2.get_connectionGenes().size());
}