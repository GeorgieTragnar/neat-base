#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <functional>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/analysis/strategies/SingleSimpleFitnessStrategy.hpp"
#include "version3/analysis/SpeciationControlUnit.hpp"
#include "version3/operator/Init.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Analysis;
using namespace Operator;

// =============================================================================
// MOCK SPECIATION CONTROL UNIT
// =============================================================================

class MockSpeciationControlUnit : public SpeciationControlUnit {
public:
    std::vector<std::shared_ptr<const Phenotype>> getChampions() const override {
        return {};
    }
    
    std::shared_ptr<const Phenotype> getBestChampion() const override {
        return nullptr;
    }
    
    std::shared_ptr<const Phenotype> getRandomChampion() const override {
        return nullptr;
    }
    
    size_t getChampionCount() const override {
        return 0;
    }
};

// =============================================================================
// TEST FITNESS RESULT IMPLEMENTATIONS
// =============================================================================

class DoubleFitnessResult : public FitnessResultInterface {
public:
    explicit DoubleFitnessResult(double value) : _value(value) {}
    
    bool isBetterThan(const FitnessResultInterface& other) const override {
        const auto* otherDouble = dynamic_cast<const DoubleFitnessResult*>(&other);
        if (!otherDouble) return false;
        return _value > otherDouble->_value;
    }
    
    bool isEqualTo(const FitnessResultInterface& other) const override {
        const auto* otherDouble = dynamic_cast<const DoubleFitnessResult*>(&other);
        if (!otherDouble) return false;
        return std::abs(_value - otherDouble->_value) < 1e-9;
    }
    
    std::unique_ptr<FitnessResultInterface> clone() const override {
        return std::make_unique<DoubleFitnessResult>(_value);
    }
    
    double getValue() const { return _value; }

private:
    double _value;
};

class IntFitnessResult : public FitnessResultInterface {
public:
    explicit IntFitnessResult(int value) : _value(value) {}
    
    bool isBetterThan(const FitnessResultInterface& other) const override {
        const auto* otherInt = dynamic_cast<const IntFitnessResult*>(&other);
        if (!otherInt) return false;
        return _value > otherInt->_value;
    }
    
    bool isEqualTo(const FitnessResultInterface& other) const override {
        const auto* otherInt = dynamic_cast<const IntFitnessResult*>(&other);
        if (!otherInt) return false;
        return _value == otherInt->_value;
    }
    
    std::unique_ptr<FitnessResultInterface> clone() const override {
        return std::make_unique<IntFitnessResult>(_value);
    }
    
    int getValue() const { return _value; }

private:
    int _value;
};

class CustomFitnessResult : public FitnessResultInterface {
public:
    explicit CustomFitnessResult(double accuracy, double speed) 
        : _accuracy(accuracy), _speed(speed) {}
    
    bool isBetterThan(const FitnessResultInterface& other) const override {
        const auto* otherCustom = dynamic_cast<const CustomFitnessResult*>(&other);
        if (!otherCustom) return false;
        // Custom fitness: weighted combination of accuracy and speed
        double myFitness = _accuracy * 0.7 + _speed * 0.3;
        double otherFitness = otherCustom->_accuracy * 0.7 + otherCustom->_speed * 0.3;
        return myFitness > otherFitness;
    }
    
    bool isEqualTo(const FitnessResultInterface& other) const override {
        const auto* otherCustom = dynamic_cast<const CustomFitnessResult*>(&other);
        if (!otherCustom) return false;
        return std::abs(_accuracy - otherCustom->_accuracy) < 1e-9 &&
               std::abs(_speed - otherCustom->_speed) < 1e-9;
    }
    
    std::unique_ptr<FitnessResultInterface> clone() const override {
        return std::make_unique<CustomFitnessResult>(_accuracy, _speed);
    }
    
    double getAccuracy() const { return _accuracy; }
    double getSpeed() const { return _speed; }

private:
    double _accuracy;
    double _speed;
};

// =============================================================================
// TEST CLASS
// =============================================================================

class SingleSimpleFitnessStrategyTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_shared<HistoryTracker>();
        mockSpeciation = std::make_shared<MockSpeciationControlUnit>();
    }

    std::shared_ptr<HistoryTracker> historyTracker;
    std::shared_ptr<MockSpeciationControlUnit> mockSpeciation;
    
    // Helper to create a simple test genome
    Genome createTestGenome() {
        std::vector<NodeGeneAttributes> inputAttrs = {{ActivationType::NONE}};
        std::vector<NodeGeneAttributes> outputAttrs = {{ActivationType::SIGMOID}};
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
        
        InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                         InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
        
        return init(historyTracker, params);
    }
    
    // Helper to create a test phenotype
    std::shared_ptr<const Phenotype> createTestPhenotype() {
        Genome genome = createTestGenome();
        genome.constructPhenotype();
        return genome.get_phenotype();
    }
};

// =============================================================================
// TEST CASES
// =============================================================================

TEST_F(SingleSimpleFitnessStrategyTest, BasicEvaluationTest) {
    // Create a simple evaluation function that returns a known value
    auto evaluationFunction = [](const Phenotype& phenotype) -> DoubleFitnessResult {
        // Simple fitness based on number of connections
        double fitness = static_cast<double>(phenotype._orderedConnections.size());
        return DoubleFitnessResult(fitness);
    };
    
    // Create strategy with the evaluation function
    SingleSimpleFitnessParams<DoubleFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<DoubleFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Evaluate using the strategy
    DoubleFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify the result matches what our evaluation function should return
    double expectedFitness = static_cast<double>(phenotype->_orderedConnections.size());
    EXPECT_EQ(result.getValue(), expectedFitness);
}

TEST_F(SingleSimpleFitnessStrategyTest, BasicEvaluationWithTestFitnessResult) {
    // Test with the common TestFitnessResult class
    auto evaluationFunction = [](const Phenotype& phenotype) -> TestFitnessResult {
        // Fitness based on number of nodes
        double fitness = static_cast<double>(phenotype._nodeGeneAttributes.size());
        return TestFitnessResult(fitness);
    };
    
    // Create strategy with TestFitnessResult
    SingleSimpleFitnessParams<TestFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<TestFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Evaluate using the strategy
    TestFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify the result
    double expectedFitness = static_cast<double>(phenotype->_nodeGeneAttributes.size());
    EXPECT_EQ(result.getFitness(), expectedFitness);
}

TEST_F(SingleSimpleFitnessStrategyTest, BasicEvaluationWithIntFitness) {
    // Test with integer-based fitness
    auto evaluationFunction = [](const Phenotype& phenotype) -> IntFitnessResult {
        // Count enabled connections
        int enabledConnections = 0;
        for (const auto& conn : phenotype._orderedConnections) {
            if (conn._connectionGeneAttribute.enabled) {
                enabledConnections++;
            }
        }
        return IntFitnessResult(enabledConnections);
    };
    
    // Create strategy with IntFitnessResult
    SingleSimpleFitnessParams<IntFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<IntFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Evaluate using the strategy
    IntFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify the result - count enabled connections manually
    int expectedEnabledConnections = 0;
    for (const auto& conn : phenotype->_orderedConnections) {
        if (conn._connectionGeneAttribute.enabled) {
            expectedEnabledConnections++;
        }
    }
    EXPECT_EQ(result.getValue(), expectedEnabledConnections);
}

TEST_F(SingleSimpleFitnessStrategyTest, BasicEvaluationWithCustomFitness) {
    // Test with complex custom fitness
    auto evaluationFunction = [](const Phenotype& phenotype) -> CustomFitnessResult {
        // Simulate complex fitness calculation
        double accuracy = static_cast<double>(phenotype._nodeGeneAttributes.size()) / 10.0;
        double speed = 1.0 / (1.0 + static_cast<double>(phenotype._orderedConnections.size()));
        return CustomFitnessResult(accuracy, speed);
    };
    
    // Create strategy with CustomFitnessResult
    SingleSimpleFitnessParams<CustomFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<CustomFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Evaluate using the strategy
    CustomFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify the result matches expected calculation
    double expectedAccuracy = static_cast<double>(phenotype->_nodeGeneAttributes.size()) / 10.0;
    double expectedSpeed = 1.0 / (1.0 + static_cast<double>(phenotype->_orderedConnections.size()));
    
    EXPECT_DOUBLE_EQ(result.getAccuracy(), expectedAccuracy);
    EXPECT_DOUBLE_EQ(result.getSpeed(), expectedSpeed);
}

TEST_F(SingleSimpleFitnessStrategyTest, EvaluationFunctionCalledTest) {
    // Test that the evaluation function is called with the correct phenotype
    bool functionWasCalled = false;
    const Phenotype* receivedPhenotype = nullptr;
    
    auto evaluationFunction = [&functionWasCalled, &receivedPhenotype](const Phenotype& phenotype) -> TestFitnessResult {
        functionWasCalled = true;
        receivedPhenotype = &phenotype;
        return TestFitnessResult(42.0);
    };
    
    // Create strategy
    SingleSimpleFitnessParams<TestFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<TestFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Verify function hasn't been called yet
    EXPECT_FALSE(functionWasCalled);
    EXPECT_EQ(receivedPhenotype, nullptr);
    
    // Evaluate using the strategy
    TestFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify the function was called with correct phenotype
    EXPECT_TRUE(functionWasCalled);
    EXPECT_EQ(receivedPhenotype, phenotype.get());
    EXPECT_EQ(result.getFitness(), 42.0);
}

TEST_F(SingleSimpleFitnessStrategyTest, ReturnValueCorrectnessTest) {
    // Test that strategy returns exactly what evaluation function returns
    struct TestData {
        DoubleFitnessResult expectedResult;
        std::string description;
    };
    
    std::vector<TestData> testCases = {
        {DoubleFitnessResult(0.0), "zero fitness"},
        {DoubleFitnessResult(1.0), "unit fitness"},
        {DoubleFitnessResult(-5.5), "negative fitness"},
        {DoubleFitnessResult(999.999), "large fitness"},
        {DoubleFitnessResult(1e-10), "very small positive fitness"}
    };
    
    for (const auto& testCase : testCases) {
        // Create evaluation function that returns specific value
        auto evaluationFunction = [&testCase](const Phenotype& phenotype) -> DoubleFitnessResult {
            return testCase.expectedResult;
        };
        
        // Create strategy
        SingleSimpleFitnessParams<DoubleFitnessResult> params{evaluationFunction};
        SingleSimpleFitnessStrategy<DoubleFitnessResult> strategy(params);
        
        // Create test phenotype
        auto phenotype = createTestPhenotype();
        
        // Evaluate using the strategy
        DoubleFitnessResult result = strategy.evaluate(*phenotype, *mockSpeciation);
        
        // Verify exact return value
        EXPECT_DOUBLE_EQ(result.getValue(), testCase.expectedResult.getValue()) 
            << "Failed for test case: " << testCase.description;
    }
}

TEST_F(SingleSimpleFitnessStrategyTest, SpeciationIgnoredTest) {
    // Test that speciation parameter doesn't affect single simple strategy results
    int evaluationCallCount = 0;
    
    auto evaluationFunction = [&evaluationCallCount](const Phenotype& phenotype) -> IntFitnessResult {
        evaluationCallCount++;
        // Return fitness based only on phenotype, ignoring any speciation context
        return IntFitnessResult(static_cast<int>(phenotype._nodeGeneAttributes.size()));
    };
    
    // Create strategy
    SingleSimpleFitnessParams<IntFitnessResult> params{evaluationFunction};
    SingleSimpleFitnessStrategy<IntFitnessResult> strategy(params);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Test with the mock speciation (should be ignored)
    IntFitnessResult result1 = strategy.evaluate(*phenotype, *mockSpeciation);
    EXPECT_EQ(evaluationCallCount, 1);
    
    // Test with a different mock speciation - result should be identical
    MockSpeciationControlUnit differentMockSpeciation;
    IntFitnessResult result2 = strategy.evaluate(*phenotype, differentMockSpeciation);
    EXPECT_EQ(evaluationCallCount, 2);
    
    // Results should be identical regardless of speciation parameter
    EXPECT_EQ(result1.getValue(), result2.getValue());
    
    // Both results should match the expected fitness (node count)
    int expectedFitness = static_cast<int>(phenotype->_nodeGeneAttributes.size());
    EXPECT_EQ(result1.getValue(), expectedFitness);
    EXPECT_EQ(result2.getValue(), expectedFitness);
}

TEST_F(SingleSimpleFitnessStrategyTest, ExceptionPropagationTest) {
    // Test that exceptions from evaluation function propagate correctly
    
    // Test case 1: Standard exception
    {
        auto evaluationFunction = [](const Phenotype& phenotype) -> TestFitnessResult {
            throw std::runtime_error("Test exception from evaluation function");
        };
        
        SingleSimpleFitnessParams<TestFitnessResult> params{evaluationFunction};
        SingleSimpleFitnessStrategy<TestFitnessResult> strategy(params);
        
        auto phenotype = createTestPhenotype();
        
        EXPECT_THROW({
            strategy.evaluate(*phenotype, *mockSpeciation);
        }, std::runtime_error);
    }
    
    // Test case 2: Custom exception
    {
        class CustomEvaluationException : public std::exception {
        public:
            const char* what() const noexcept override {
                return "Custom evaluation exception";
            }
        };
        
        auto evaluationFunction = [](const Phenotype& phenotype) -> DoubleFitnessResult {
            throw CustomEvaluationException();
        };
        
        SingleSimpleFitnessParams<DoubleFitnessResult> params{evaluationFunction};
        SingleSimpleFitnessStrategy<DoubleFitnessResult> strategy(params);
        
        auto phenotype = createTestPhenotype();
        
        EXPECT_THROW({
            strategy.evaluate(*phenotype, *mockSpeciation);
        }, CustomEvaluationException);
    }
    
    // Test case 3: Verify strategy doesn't catch or modify exceptions
    {
        bool exceptionWasThrown = false;
        auto evaluationFunction = [&exceptionWasThrown](const Phenotype& phenotype) -> IntFitnessResult {
            exceptionWasThrown = true;
            throw std::invalid_argument("Invalid phenotype data");
        };
        
        SingleSimpleFitnessParams<IntFitnessResult> params{evaluationFunction};
        SingleSimpleFitnessStrategy<IntFitnessResult> strategy(params);
        
        auto phenotype = createTestPhenotype();
        
        try {
            strategy.evaluate(*phenotype, *mockSpeciation);
            FAIL() << "Expected exception was not thrown";
        } catch (const std::invalid_argument& e) {
            EXPECT_TRUE(exceptionWasThrown);
            EXPECT_STREQ(e.what(), "Invalid phenotype data");
        } catch (...) {
            FAIL() << "Unexpected exception type was thrown";
        }
    }
}

TEST_F(SingleSimpleFitnessStrategyTest, TemplateInstantiationTest) {
    // Test that multiple template instantiations of SingleSimpleFitnessStrategy can coexist
    // This verifies the template system works correctly with different fitness result types
    
    // Create evaluation functions for different fitness result types
    auto doubleEvaluationFunction = [](const Phenotype& phenotype) -> DoubleFitnessResult {
        return DoubleFitnessResult(static_cast<double>(phenotype._nodeGeneAttributes.size()));
    };
    
    auto intEvaluationFunction = [](const Phenotype& phenotype) -> IntFitnessResult {
        return IntFitnessResult(static_cast<int>(phenotype._orderedConnections.size()));
    };
    
    auto testEvaluationFunction = [](const Phenotype& phenotype) -> TestFitnessResult {
        return TestFitnessResult(42.0);
    };
    
    auto customEvaluationFunction = [](const Phenotype& phenotype) -> CustomFitnessResult {
        return CustomFitnessResult(0.8, 0.6);
    };
    
    // Create multiple strategy instances with different template types
    SingleSimpleFitnessParams<DoubleFitnessResult> doubleParams{doubleEvaluationFunction};
    SingleSimpleFitnessStrategy<DoubleFitnessResult> doubleStrategy(doubleParams);
    
    SingleSimpleFitnessParams<IntFitnessResult> intParams{intEvaluationFunction};
    SingleSimpleFitnessStrategy<IntFitnessResult> intStrategy(intParams);
    
    SingleSimpleFitnessParams<TestFitnessResult> testParams{testEvaluationFunction};
    SingleSimpleFitnessStrategy<TestFitnessResult> testStrategy(testParams);
    
    SingleSimpleFitnessParams<CustomFitnessResult> customParams{customEvaluationFunction};
    SingleSimpleFitnessStrategy<CustomFitnessResult> customStrategy(customParams);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Evaluate using all different strategy types
    DoubleFitnessResult doubleResult = doubleStrategy.evaluate(*phenotype, *mockSpeciation);
    IntFitnessResult intResult = intStrategy.evaluate(*phenotype, *mockSpeciation);
    TestFitnessResult testResult = testStrategy.evaluate(*phenotype, *mockSpeciation);
    CustomFitnessResult customResult = customStrategy.evaluate(*phenotype, *mockSpeciation);
    
    // Verify each strategy returns the expected type and value
    EXPECT_EQ(doubleResult.getValue(), static_cast<double>(phenotype->_nodeGeneAttributes.size()));
    EXPECT_EQ(intResult.getValue(), static_cast<int>(phenotype->_orderedConnections.size()));
    EXPECT_EQ(testResult.getFitness(), 42.0);
    EXPECT_DOUBLE_EQ(customResult.getAccuracy(), 0.8);
    EXPECT_DOUBLE_EQ(customResult.getSpeed(), 0.6);
    
    // Verify that all template instantiations can coexist in the same compilation unit
    // by storing them in containers and using them together
    std::vector<std::function<void()>> evaluators = {
        [&]() { doubleStrategy.evaluate(*phenotype, *mockSpeciation); },
        [&]() { intStrategy.evaluate(*phenotype, *mockSpeciation); },
        [&]() { testStrategy.evaluate(*phenotype, *mockSpeciation); },
        [&]() { customStrategy.evaluate(*phenotype, *mockSpeciation); }
    };
    
    // Execute all evaluators to ensure no compilation or linking issues
    for (auto& evaluator : evaluators) {
        EXPECT_NO_THROW(evaluator());
    }
}

TEST_F(SingleSimpleFitnessStrategyTest, DefaultConstructedParamsTest) {
    // Test behavior with default-constructed parameters (null evaluation function)
    SingleSimpleFitnessParams<TestFitnessResult> defaultParams{};
    
    // The default-constructed std::function should be empty/null
    EXPECT_FALSE(defaultParams.evaluationFunction) << "Default-constructed evaluation function should be null";
    
    // Creating strategy with null evaluation function
    SingleSimpleFitnessStrategy<TestFitnessResult> strategy(defaultParams);
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Calling evaluate with null evaluation function should throw or have undefined behavior
    // std::function throws std::bad_function_call when called while empty
    EXPECT_THROW({
        strategy.evaluate(*phenotype, *mockSpeciation);
    }, std::bad_function_call) << "Evaluating with null evaluation function should throw std::bad_function_call";
}

TEST_F(SingleSimpleFitnessStrategyTest, PolymorphicUsageTest) {
    // Test using strategy through base class pointer to verify virtual function override
    auto evaluationFunction = [](const Phenotype& phenotype) -> TestFitnessResult {
        return TestFitnessResult(123.45);
    };
    
    SingleSimpleFitnessParams<TestFitnessResult> params{evaluationFunction};
    
    // Create strategy through base class pointer
    std::unique_ptr<FitnessStrategy<TestFitnessResult>> strategy = 
        std::make_unique<SingleSimpleFitnessStrategy<TestFitnessResult>>(params);
    
    // Verify the strategy object was created successfully
    ASSERT_NE(strategy, nullptr) << "Strategy should be successfully created through base class pointer";
    
    // Create test phenotype
    auto phenotype = createTestPhenotype();
    
    // Call evaluate through base class interface
    TestFitnessResult result = strategy->evaluate(*phenotype, *mockSpeciation);
    
    // Verify the virtual function call worked correctly
    EXPECT_EQ(result.getFitness(), 123.45) << "Polymorphic call should return expected fitness value";
    
    // Test with different strategy instance to verify polymorphic behavior
    auto differentEvaluationFunction = [](const Phenotype& phenotype) -> TestFitnessResult {
        return TestFitnessResult(-99.0);
    };
    
    SingleSimpleFitnessParams<TestFitnessResult> differentParams{differentEvaluationFunction};
    std::unique_ptr<FitnessStrategy<TestFitnessResult>> differentStrategy = 
        std::make_unique<SingleSimpleFitnessStrategy<TestFitnessResult>>(differentParams);
    
    TestFitnessResult differentResult = differentStrategy->evaluate(*phenotype, *mockSpeciation);
    EXPECT_EQ(differentResult.getFitness(), -99.0) << "Different polymorphic instance should return different result";
    
    // Verify strategies are independent
    TestFitnessResult originalResult = strategy->evaluate(*phenotype, *mockSpeciation);
    EXPECT_EQ(originalResult.getFitness(), 123.45) << "Original strategy should be unaffected by other instances";
}
