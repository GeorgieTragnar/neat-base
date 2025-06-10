# NEAT Testing Framework

This directory contains the comprehensive testing suite for the NEAT (NeuroEvolution of Augmenting Topologies) implementation, built on Google Test with custom infrastructure for automated testing of genetic algorithm components.

## Architecture Overview

### Test Organization
```
tests/
├── main.cpp                    # Test runner entry point
├── test_common.h              # Global testing utilities and access overrides
├── test_utilities.h/cpp       # Shared testing infrastructure
├── unit/                      # Unit tests organized by component
│   ├── ConnectionGene/        # Connection gene tests with code generation
│   ├── NodeGene/             # Node gene tests with code generation  
│   ├── Genome/               # Core genome data structure tests
│   ├── Init/                 # Initialization operator tests
│   ├── ConnectionMutation/   # Connection mutation operator tests
│   ├── NodeMutation/         # Node mutation operator tests
│   ├── WeightMutation/       # Weight mutation operator tests
│   ├── ConnectionReactivation/ # Connection reactivation tests
│   ├── Crossover/            # Crossover operator tests
│   ├── CycleDetection/       # Cycle detection tests
│   ├── HistoryTracker/       # History tracking tests
│   ├── SingleSimpleFitnessStrategy/ # Fitness strategy tests
│   └── Visualization/        # Visualization tests
└── integration/              # Integration tests (future expansion)
```

## Build System

### CMake Configuration
The testing framework uses a sophisticated CMake build system with:

- **Automatic Subdirectory Discovery**: Processes subdirectories with specialized build requirements
- **File Collection**: Aggregates test files and include directories using parent scope variables
- **Unified Executable**: Creates single `neat_tests` executable containing all tests
- **Code Generation Integration**: Supports automatic code generation for attribute testing

### Build Commands
```bash
# From build directory
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .

# Run tests
ctest
# or
./neat_tests
```

### Specialized Build Patterns

**Components with Code Generation** (NodeGene, ConnectionGene):
```cmake
include(generate_*_registration.cmake)
generate_*_registration(source_file output_file)
add_custom_command() # Regenerate when headers change
```

**Simple Components** (Genome, operators):
```cmake
set(COMPONENT_TEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/Component.cpp PARENT_SCOPE)
```

## Testing Infrastructure

### Core Components

**test_common.h**
- Redefines `protected` and `private` as `public` when `TESTING` is defined
- Enables white-box testing of internal class members
- Must be included before project headers

**test_utilities.h/cpp**
- Generic testing framework with runtime type introspection
- Comprehensive value generation for any type
- Type-safe comparison and modification functions
- Automatic attribute combination testing

### Key Utilities

**Generic Value Generation:**
```cpp
// Automatically generates test values for any type
auto values = generate_runtime_test_values<T>();

// Creates all combinations of attribute values  
auto combinations = generate_runtime_attribute_combinations<Attributes, Registry>();
```

**Type-Safe Operations:**
```cpp
// Generic comparison with floating-point tolerance
bool equal = runtime_attributes_equal<Attributes, Registry>(a, b);

// Safe attribute modification for testing
runtime_modify_attributes_systematically<Attributes, Registry>(attrs, field_index);
```

## Code Generation System

### Automatic Field Registration
- **Parser**: Regex-based parsing of struct definitions from header files
- **Generated Files**: `.inc` files with field registration code
- **Build Integration**: Custom commands ensure regeneration when headers change
- **Runtime Registry**: Dynamic field descriptors with type information and operations

**Field Descriptor Structure:**
```cpp
struct FieldDescriptor {
    std::string name;
    std::type_index type;
    std::function<std::vector<std::any>()> value_generator;
    std::function<void(Attributes&, const std::any&)> setter;
    std::function<std::any(const Attributes&)> getter;
    bool is_primitive;
};
```

## Testing Patterns

### 1. Construction Testing
- **Exhaustive Combinations**: Tests all combinations of valid parameters
- **Validation Checking**: Verifies death tests for invalid parameters
- **Immutability Contracts**: Ensures certain fields remain immutable after construction

Example:
```cpp
TEST_F(ComponentTest, ConstructorCombinations) {
    auto combinations = generate_runtime_attribute_combinations<Attributes, Registry>();
    for (const auto& combo : combinations) {
        // Test construction with each combination
        Component component(combo);
        EXPECT_TRUE(component.isValid());
    }
}
```

### 2. Copy Semantics Testing
- **Exact Replication**: Copy constructors produce identical objects
- **Independence**: Modifications to copies don't affect originals
- **Reference Remapping**: Complex object graphs maintain correct references
- **Self-Assignment**: Assignment operators handle self-assignment correctly

Example:
```cpp
TEST_F(GenomeTest, CopyConstructorIndependence) {
    Genome original = createTestGenome();
    Genome copy(original);
    
    // Modify copy
    copy.mutate();
    
    // Original should be unchanged
    EXPECT_NE(original, copy);
}
```

### 3. Runtime Extensibility
- **Automatic Adaptation**: Tests automatically adapt when struct fields change
- **Field Validation**: Runtime checks ensure all fields are primitive types
- **Size Monitoring**: Tracks struct sizes to detect unexpected bloat
- **Combination Generation**: Creates comprehensive test scenarios

### 4. Operator Testing
- **Parameter Validation**: Tests all parameter combinations
- **Reproducibility**: Verifies deterministic behavior with fixed seeds
- **Edge Cases**: Tests boundary conditions and invalid inputs
- **Mutation Coverage**: Ensures all possible mutations are tested

Example:
```cpp
TEST_F(MutationTest, ReproducibilityWithFixedSeed) {
    auto params = MutationParams{/* ... */};
    
    // Same seed should produce same results
    auto result1 = mutation_operator(genome, params, 12345);
    auto result2 = mutation_operator(genome, params, 12345);
    
    EXPECT_EQ(result1, result2);
}
```

## Advanced Features

### Death Testing
Extensive use of `EXPECT_DEATH()` for validation testing:
```cpp
TEST_F(ComponentTest, InvalidParametersDeath) {
    EXPECT_DEATH(Component(invalid_params), "Assertion failed");
}
```

### Mock Implementations
Test-specific implementations for complex dependencies:
```cpp
class TestFitnessResult : public FitnessResult {
    // Simplified implementation for testing
};
```

### Performance Testing
- Test count limiting to avoid excessive combinations
- Efficient type dispatch systems
- Minimal overhead generic testing infrastructure

## Usage Guidelines

### Running Tests
```bash
# Run all tests
./neat_tests

# Run specific test suite
./neat_tests --gtest_filter="GenomeTest.*"

# Run with verbose output
./neat_tests --gtest_filter="*" --gtest_print_time=1
```

### Adding New Tests

**For Simple Components:**
1. Create `tests/unit/ComponentName/ComponentName.cpp`
2. Add CMakeLists.txt if specialized build requirements exist
3. Follow existing patterns for construction, copy semantics, and attribute testing

**For Components with Attributes:**
1. Ensure struct fields are properly defined in header
2. Add code generation support if needed
3. Use runtime attribute testing utilities
4. Test all attribute combinations systematically

**Test Structure Template:**
```cpp
#include "test_common.h"
#include "test_utilities.h"
#include "version3/path/to/Component.hpp"

class ComponentTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup test fixtures
    }
    
    // Helper methods
};

TEST_F(ComponentTest, Construction) {
    // Test construction patterns
}

TEST_F(ComponentTest, CopySemantics) {
    // Test copy constructor and assignment
}

TEST_F(ComponentTest, AttributeModification) {
    // Test attribute modification if applicable
}
```

## Dependencies

- **Google Test**: Primary testing framework
- **nlohmann_json**: JSON handling for test data
- **spdlog**: Logging within tests
- **Standard C++20**: Template metaprogramming and concepts

## Best Practices

1. **Automatic Adaptation**: Write tests that adapt to code changes automatically
2. **Exhaustive Coverage**: Test all valid parameter combinations
3. **Type Safety**: Use generic testing infrastructure that maintains type safety
4. **Clear Naming**: Use descriptive test names that explain what is being tested
5. **Focused Tests**: Each test should verify one specific behavior
6. **Proper Fixtures**: Use test fixtures for complex setup and teardown
7. **Documentation**: Comment complex test logic and edge cases

## Contributing

When adding new tests:
1. Follow existing patterns and conventions
2. Use the generic testing infrastructure when possible
3. Add appropriate death tests for invalid inputs
4. Verify tests pass with both debug and release builds
5. Update this README if introducing new testing patterns

The testing framework is designed to be maintainable and extensible, automatically adapting to changes in the codebase while providing comprehensive coverage of the NEAT implementation.