#include "tests/test_utilities.h"
#include "version3/data/NodeGene.hpp"
#include <stdexcept>
#include <string>

// =============================================================================
// TYPE DISPATCH FUNCTIONS
// =============================================================================

bool dispatch_comparison(const std::type_index& type_idx, const std::any& a, const std::any& b) {
    // Handle common types efficiently
    if (type_idx == std::type_index(typeid(bool))) {
        return compare_values<bool>(a, b);
    }
    else if (type_idx == std::type_index(typeid(int))) {
        return compare_values<int>(a, b);
    }
    else if (type_idx == std::type_index(typeid(float))) {
        return compare_values<float>(a, b);
    }
    else if (type_idx == std::type_index(typeid(double))) {
        return compare_values<double>(a, b);
    }
    // Integer types
    else if (type_idx == std::type_index(typeid(int8_t))) {
        return compare_values<int8_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(uint8_t))) {
        return compare_values<uint8_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(int16_t))) {
        return compare_values<int16_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(uint16_t))) {
        return compare_values<uint16_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(int32_t))) {
        return compare_values<int32_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(uint32_t))) {
        return compare_values<uint32_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(int64_t))) {
        return compare_values<int64_t>(a, b);
    }
    else if (type_idx == std::type_index(typeid(uint64_t))) {
        return compare_values<uint64_t>(a, b);
    }
    // Additional integer types
    else if (type_idx == std::type_index(typeid(char))) {
        return compare_values<char>(a, b);
    }
    else if (type_idx == std::type_index(typeid(signed char))) {
        return compare_values<signed char>(a, b);
    }
    else if (type_idx == std::type_index(typeid(unsigned char))) {
        return compare_values<unsigned char>(a, b);
    }
    else if (type_idx == std::type_index(typeid(short))) {
        return compare_values<short>(a, b);
    }
    else if (type_idx == std::type_index(typeid(unsigned short))) {
        return compare_values<unsigned short>(a, b);
    }
    else if (type_idx == std::type_index(typeid(long))) {
        return compare_values<long>(a, b);
    }
    else if (type_idx == std::type_index(typeid(unsigned long))) {
        return compare_values<unsigned long>(a, b);
    }
    else if (type_idx == std::type_index(typeid(long long))) {
        return compare_values<long long>(a, b);
    }
    else if (type_idx == std::type_index(typeid(unsigned long long))) {
        return compare_values<unsigned long long>(a, b);
    }
    // Floating point types
    else if (type_idx == std::type_index(typeid(long double))) {
        return compare_values<long double>(a, b);
    }
    // Enum types (add more as needed)
    else if (type_idx == std::type_index(typeid(ActivationType))) {
        return compare_values<ActivationType>(a, b);
    }
    else if (type_idx == std::type_index(typeid(NodeType))) {
        return compare_values<NodeType>(a, b);
    }
    else {
        throw std::runtime_error("Unsupported field type for comparison: " + std::string(type_idx.name()) + 
                                 ". Add this type to dispatch_comparison function.");
    }
}

bool dispatch_modification(const std::type_index& type_idx, std::any& current_value, const std::vector<std::any>& test_values) {
    // Handle common types efficiently
    if (type_idx == std::type_index(typeid(bool))) {
        return modify_value<bool>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(int))) {
        return modify_value<int>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(float))) {
        return modify_value<float>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(double))) {
        return modify_value<double>(current_value, test_values);
    }
    // Integer types
    else if (type_idx == std::type_index(typeid(int8_t))) {
        return modify_value<int8_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(uint8_t))) {
        return modify_value<uint8_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(int16_t))) {
        return modify_value<int16_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(uint16_t))) {
        return modify_value<uint16_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(int32_t))) {
        return modify_value<int32_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(uint32_t))) {
        return modify_value<uint32_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(int64_t))) {
        return modify_value<int64_t>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(uint64_t))) {
        return modify_value<uint64_t>(current_value, test_values);
    }
    // Additional integer types
    else if (type_idx == std::type_index(typeid(char))) {
        return modify_value<char>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(signed char))) {
        return modify_value<signed char>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(unsigned char))) {
        return modify_value<unsigned char>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(short))) {
        return modify_value<short>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(unsigned short))) {
        return modify_value<unsigned short>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(long))) {
        return modify_value<long>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(unsigned long))) {
        return modify_value<unsigned long>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(long long))) {
        return modify_value<long long>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(unsigned long long))) {
        return modify_value<unsigned long long>(current_value, test_values);
    }
    // Floating point types
    else if (type_idx == std::type_index(typeid(long double))) {
        return modify_value<long double>(current_value, test_values);
    }
    // Enum types (add more as needed)
    else if (type_idx == std::type_index(typeid(ActivationType))) {
        return modify_value<ActivationType>(current_value, test_values);
    }
    else if (type_idx == std::type_index(typeid(NodeType))) {
        return modify_value<NodeType>(current_value, test_values);
    }
    else {
        throw std::runtime_error("Unsupported field type for modification: " + std::string(type_idx.name()) + 
                                 ". Add this type to dispatch_modification function.");
    }
}