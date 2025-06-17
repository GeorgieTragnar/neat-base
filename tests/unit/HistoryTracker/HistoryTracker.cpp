#include <gtest/gtest.h>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/data/HistoryTracker.hpp"
#include "version3/operator/Init.hpp"

class HistoryTrackerTest : public ::testing::Test {
protected:
    HistoryTracker tracker;
};

TEST_F(HistoryTrackerTest, NodeIDsIncrement) {
    uint32_t first = tracker.get_input(0);
    uint32_t second = tracker.get_input(1);
    EXPECT_EQ(first, 1u);
    EXPECT_EQ(second, 2u);
}

TEST_F(HistoryTrackerTest, ConnectionIDsIncrement) {
    uint32_t first = tracker.get_connection(1, 2);
    uint32_t second = tracker.get_connection(3, 4);
    EXPECT_EQ(first, 1u);
    EXPECT_EQ(second, 2u);
}

TEST_F(HistoryTrackerTest, InputConsistency) {
    uint32_t first = tracker.get_input(0);
    uint32_t second = tracker.get_input(0);
    EXPECT_EQ(first, second);
}

TEST_F(HistoryTrackerTest, InputUniqueness) {
    uint32_t input0 = tracker.get_input(0);
    uint32_t input1 = tracker.get_input(1);
    EXPECT_NE(input0, input1);
}

TEST_F(HistoryTrackerTest, OutputConsistency) {
    uint32_t first = tracker.get_output(0);
    uint32_t second = tracker.get_output(0);
    EXPECT_EQ(first, second);
}

TEST_F(HistoryTrackerTest, OutputUniqueness) {
    uint32_t output0 = tracker.get_output(0);
    uint32_t output1 = tracker.get_output(1);
    EXPECT_NE(output0, output1);
}

TEST_F(HistoryTrackerTest, SparseInputAccess) {
    uint32_t input5 = tracker.get_input(5);
    uint32_t input2 = tracker.get_input(2);
    uint32_t input5_again = tracker.get_input(5);
    uint32_t input2_again = tracker.get_input(2);
    
    EXPECT_EQ(input5, input5_again);
    EXPECT_EQ(input2, input2_again);
    EXPECT_NE(input5, input2);
}

TEST_F(HistoryTrackerTest, SparseOutputAccess) {
    uint32_t output3 = tracker.get_output(3);
    uint32_t output1 = tracker.get_output(1);
    uint32_t output3_again = tracker.get_output(3);
    uint32_t output1_again = tracker.get_output(1);
    
    EXPECT_EQ(output3, output3_again);
    EXPECT_EQ(output1, output1_again);
    EXPECT_NE(output3, output1);
}

TEST_F(HistoryTrackerTest, BiasConsistency) {
    uint32_t first = tracker.get_bias();
    uint32_t second = tracker.get_bias();
    EXPECT_EQ(first, second);
}

TEST_F(HistoryTrackerTest, BiasUniqueness) {
    uint32_t bias = tracker.get_bias();
    uint32_t input = tracker.get_input(0);
    uint32_t output = tracker.get_output(0);
    
    EXPECT_NE(bias, input);
    EXPECT_NE(bias, output);
    EXPECT_NE(input, output);
}

TEST_F(HistoryTrackerTest, ConnectionConsistency) {
    uint32_t first = tracker.get_connection(10, 20);
    uint32_t second = tracker.get_connection(10, 20);
    EXPECT_EQ(first, second);
}

TEST_F(HistoryTrackerTest, ConnectionOrderMatters) {
    uint32_t ab = tracker.get_connection(10, 20);
    uint32_t ba = tracker.get_connection(20, 10);
    EXPECT_NE(ab, ba);
}

TEST_F(HistoryTrackerTest, ConnectionUniqueness) {
    uint32_t ab = tracker.get_connection(10, 20);
    uint32_t ac = tracker.get_connection(10, 30);
    EXPECT_NE(ab, ac);
}

TEST_F(HistoryTrackerTest, SplitNodeConsistency) {
    uint32_t connID = tracker.get_connection(1, 2);
    uint32_t first = tracker.get_splitNode(connID);
    uint32_t second = tracker.get_splitNode(connID);
    EXPECT_EQ(first, second);
}

TEST_F(HistoryTrackerTest, SplitBranchUniqueness) {
    uint32_t connID = tracker.get_connection(1, 2);
    uint32_t branch1 = tracker.create_splitBranch(connID);
    uint32_t branch2 = tracker.create_splitBranch(connID);
    EXPECT_NE(branch1, branch2);
}

TEST_F(HistoryTrackerTest, GetAllSplitNodesBasic) {
    uint32_t connID = tracker.get_connection(1, 2);
    uint32_t primaryNode = tracker.get_splitNode(connID);
    
    std::vector<uint32_t> allNodes = tracker.get_allSplitNodes(connID);
    EXPECT_EQ(allNodes.size(), 1u);
    EXPECT_EQ(allNodes[0], primaryNode);
}

TEST_F(HistoryTrackerTest, ReactivationScenario) {
    uint32_t connID = tracker.get_connection(1, 2);
    uint32_t primaryNode = tracker.get_splitNode(connID);
    uint32_t branch1 = tracker.create_splitBranch(connID);
    uint32_t branch2 = tracker.create_splitBranch(connID);
    
    std::vector<uint32_t> allNodes = tracker.get_allSplitNodes(connID);
    EXPECT_EQ(allNodes.size(), 3u);
    EXPECT_EQ(allNodes[0], primaryNode);
    EXPECT_EQ(allNodes[1], branch1);
    EXPECT_EQ(allNodes[2], branch2);
}

TEST_F(HistoryTrackerTest, SplitBeforePrimary) {
    uint32_t connID = tracker.get_connection(1, 2);
    uint32_t branch1 = tracker.create_splitBranch(connID);
    uint32_t branch2 = tracker.create_splitBranch(connID);
    uint32_t primaryNode = tracker.get_splitNode(connID);
    
    std::vector<uint32_t> allNodes = tracker.get_allSplitNodes(connID);
    EXPECT_EQ(allNodes.size(), 3u);
    EXPECT_EQ(allNodes[0], primaryNode);
    EXPECT_EQ(allNodes[1], branch1);
    EXPECT_EQ(allNodes[2], branch2);
}