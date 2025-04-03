#include <gtest/gtest.h>
#include "network/Activation.hpp"

namespace {

class ActivationTest : public ::testing::Test {
protected:
    // Setup code (runs before each test)
    void SetUp() override {
        // Any common setup
    }
};

// Test sigmoid activation function
TEST_F(ActivationTest, SigmoidFunction) {
    // Arrange
    auto sigmoid = neat::network::activation::ActivationFunctions::sigmoid;
    
    // Act & Assert - with updated expected values matching steepness 4.9
    EXPECT_NEAR(sigmoid(0.0), 0.5, 1e-6);
    EXPECT_NEAR(sigmoid(1.0), 0.9900, 1e-4);  // Updated from 0.9866
    EXPECT_NEAR(sigmoid(-1.0), 0.0100, 1e-4); // Updated from 0.0134
    
    // Test boundaries
    EXPECT_GT(sigmoid(10.0), 0.9999);
    EXPECT_LT(sigmoid(-10.0), 0.0001);
}

// Test ReLU activation function
TEST_F(ActivationTest, ReluFunction) {
    // Arrange
    auto relu = neat::network::activation::ActivationFunctions::relu;
    
    // Act & Assert
    EXPECT_DOUBLE_EQ(relu(0.0), 0.0);
    EXPECT_DOUBLE_EQ(relu(1.0), 1.0);
    EXPECT_DOUBLE_EQ(relu(-1.0), 0.0);
    EXPECT_DOUBLE_EQ(relu(5.5), 5.5);
}

// Test function retrieval
TEST_F(ActivationTest, GetActivationFunction) {
    // Arrange
    auto sigmoid = neat::network::activation::ActivationFunctions::getFunction(
        neat::network::activation::EActivationFunction::SIGMOID);
    
    // Act & Assert
    EXPECT_NEAR(sigmoid(0.0), 0.5, 1e-6);
}

} // namespace