#pragma once

#include <vector>
#include <any>
#include <typeindex>
#include <functional>
#include <limits>
#include <type_traits>
#include <memory>
#include <cmath>

#include "version3/analysis/FitnessResult.hpp"

// =============================================================================
// TEST FITNESS RESULT IMPLEMENTATION
// =============================================================================

// Test-specific fitness result implementation for use across test files
class TestFitnessResult : public Analysis::FitnessResultInterface {
public:
    explicit TestFitnessResult(double fitness) : _fitness(fitness) {}
    
    bool isBetterThan(const Analysis::FitnessResultInterface& other) const override {
        const auto* otherTest = dynamic_cast<const TestFitnessResult*>(&other);
        if (!otherTest) return false;
        return _fitness > otherTest->_fitness;
    }
    
    bool isEqualTo(const Analysis::FitnessResultInterface& other) const override {
        const auto* otherTest = dynamic_cast<const TestFitnessResult*>(&other);
        if (!otherTest) return false;
        return std::abs(_fitness - otherTest->_fitness) < 1e-9;
    }
    
    std::unique_ptr<Analysis::FitnessResultInterface> clone() const override {
        return std::make_unique<TestFitnessResult>(_fitness);
    }
    
    double getFitness() const { return _fitness; }
    
    // Comparison operator needed for std::multimap ordering
    bool operator<(const TestFitnessResult& other) const {
        return _fitness < other._fitness;
    }

private:
    double _fitness;
};

// =============================================================================
// PRIMITIVE TYPE DETECTION
// =============================================================================

template<typename T>
bool is_primitive_type() {
    return std::is_arithmetic_v<T> || std::is_enum_v<T>;
}

// =============================================================================
// RUNTIME TEST VALUE GENERATION
// =============================================================================

template<typename T>
std::vector<std::any> generate_runtime_test_values() {
    std::vector<std::any> values;
    
    if constexpr (std::is_same_v<T, bool>) {
        values.push_back(std::any(false));
        values.push_back(std::any(true));
    }
    else if constexpr (std::is_integral_v<T> && !std::is_enum_v<T>) {
        values.push_back(std::any(T{0}));
        values.push_back(std::any(T{1}));
        values.push_back(std::any(std::numeric_limits<T>::min()));
        values.push_back(std::any(std::numeric_limits<T>::max()));
        values.push_back(std::any(T{42}));
    }
    else if constexpr (std::is_floating_point_v<T>) {
        values.push_back(std::any(T{0.0}));
        values.push_back(std::any(T{1.0}));
        values.push_back(std::any(T{-1.0}));
        values.push_back(std::any(std::numeric_limits<T>::min()));
        values.push_back(std::any(std::numeric_limits<T>::max()));
        values.push_back(std::any(T{3.14159}));
        values.push_back(std::any(T{-2.5}));
    }
    else if constexpr (std::is_enum_v<T>) {
        using UnderlyingType = std::underlying_type_t<T>;
        for (UnderlyingType i = 0; i < 16; ++i) {
            values.push_back(std::any(static_cast<T>(i)));
        }
    }
    
    return values;
}

// =============================================================================
// GENERIC COMPARISON AND MODIFICATION FUNCTIONS
// =============================================================================

// Generic comparison for any two values of the same type
template<typename T>
bool compare_values(const std::any& a, const std::any& b) {
    try {
        T val_a = std::any_cast<T>(a);
        T val_b = std::any_cast<T>(b);
        
        if constexpr (std::is_floating_point_v<T>) {
            return std::abs(val_a - val_b) < std::numeric_limits<T>::epsilon();
        } else {
            return val_a == val_b;
        }
    } catch (const std::bad_any_cast&) {
        return false;
    }
}

// Generic value modification for any type
template<typename T>
bool modify_value(std::any& current_value, const std::vector<std::any>& test_values) {
    try {
        T current = std::any_cast<T>(current_value);
        
        for (const auto& test_value : test_values) {
            try {
                T test = std::any_cast<T>(test_value);
                if constexpr (std::is_floating_point_v<T>) {
                    if (std::abs(current - test) >= std::numeric_limits<T>::epsilon()) {
                        current_value = test_value;
                        return true;
                    }
                } else {
                    if (current != test) {
                        current_value = test_value;
                        return true;
                    }
                }
            } catch (const std::bad_any_cast&) {
                continue;
            }
        }
        return false;
    } catch (const std::bad_any_cast&) {
        return false;
    }
}

// =============================================================================
// TYPE DISPATCH FUNCTIONS (DECLARED HERE, DEFINED IN .cpp)
// =============================================================================

// Dispatch comparison based on type category
bool dispatch_comparison(const std::type_index& type_idx, const std::any& a, const std::any& b);

// Dispatch modification based on type category  
bool dispatch_modification(const std::type_index& type_idx, std::any& current_value, const std::vector<std::any>& test_values);

// =============================================================================
// GENERIC ATTRIBUTE COMBINATION GENERATION
// =============================================================================

template<typename AttributeType, typename RegistryType>
std::vector<AttributeType> generate_runtime_attribute_combinations() {
    const auto& fields = RegistryType::get_fields();
    std::vector<AttributeType> combinations;
    
    if (fields.empty()) {
        combinations.push_back(AttributeType{});
        return combinations;
    }
    
    // Get test values for each field
    std::vector<std::vector<std::any>> field_values;
    for (const auto& field : fields) {
        field_values.push_back(field.value_generator());
    }
    
    // Generate all combinations using recursive approach
    std::function<void(size_t, AttributeType&)> generate_combinations = 
        [&](size_t field_index, AttributeType& current_attrs) {
            if (field_index >= fields.size()) {
                combinations.push_back(current_attrs);
                return;
            }
            
            for (const auto& value : field_values[field_index]) {
                AttributeType next_attrs = current_attrs;
                fields[field_index].setter(next_attrs, value);
                generate_combinations(field_index + 1, next_attrs);
            }
        };
    
    AttributeType base_attrs{};
    generate_combinations(0, base_attrs);
    
    return combinations;
}

// =============================================================================
// GENERIC ATTRIBUTE OPERATIONS
// =============================================================================

template<typename AttributeType, typename RegistryType>
bool runtime_attributes_equal(const AttributeType& a, const AttributeType& b) {
    const auto& fields = RegistryType::get_fields();
    
    for (const auto& field : fields) {
        std::any value_a = field.getter(a);
        std::any value_b = field.getter(b);
        
        if (!dispatch_comparison(field.type, value_a, value_b)) {
            return false;
        }
    }
    
    return true;
}

template<typename AttributeType, typename RegistryType>
void runtime_modify_attributes_systematically(AttributeType& attrs, size_t field_index) {
    const auto& fields = RegistryType::get_fields();
    
    if (field_index >= fields.size()) {
        return;
    }
    
    const auto& field = fields[field_index];
    auto test_values = field.value_generator();
    
    if (test_values.empty()) {
        return;
    }
    
    std::any current_value = field.getter(attrs);
    
    if (dispatch_modification(field.type, current_value, test_values)) {
        field.setter(attrs, current_value);
    }
}