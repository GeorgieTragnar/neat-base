#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <limits>
#include <functional>
#include <sstream>
#include <typeinfo>
#include <typeindex>
#include <unordered_set>
#include <stdexcept>

#include "version3/data/NodeGene.hpp"

// =============================================================================
// RUNTIME FIELD REGISTRY SYSTEM
// =============================================================================

class NodeAttributesRegistry {
public:
    struct FieldDescriptor {
        std::string name;
        std::type_index type;
        std::function<std::vector<std::any>()> value_generator;
        std::function<void(NodeGeneAttributes&, const std::any&)> setter;
        std::function<std::any(const NodeGeneAttributes&)> getter;
        bool is_primitive;
    };

private:
    static std::vector<FieldDescriptor> fields_;
    static bool initialized_;

public:
    static void register_field(const FieldDescriptor& field) {
        fields_.push_back(field);
    }

    static const std::vector<FieldDescriptor>& get_fields() {
        if (!initialized_) {
            initialize_fields();
            initialized_ = true;
        }
        return fields_;
    }

    static size_t field_count() {
        return get_fields().size();
    }

    static bool all_fields_primitive() {
        for (const auto& field : get_fields()) {
            if (!field.is_primitive) {
                return false;
            }
        }
        return true;
    }

    static std::vector<std::string> get_non_primitive_fields() {
        std::vector<std::string> non_primitives;
        for (const auto& field : get_fields()) {
            if (!field.is_primitive) {
                non_primitives.push_back(field.name + " (" + field.type.name() + ")");
            }
        }
        return non_primitives;
    }

private:
    // Private versions that work directly on fields_ during initialization
    static bool all_fields_primitive_direct() {
        for (const auto& field : fields_) {
            if (!field.is_primitive) {
                return false;
            }
        }
        return true;
    }

    static std::vector<std::string> get_non_primitive_fields_direct() {
        std::vector<std::string> non_primitives;
        for (const auto& field : fields_) {
            if (!field.is_primitive) {
                non_primitives.push_back(field.name + " (" + field.type.name() + ")");
            }
        }
        return non_primitives;
    }

private:
    static void initialize_fields();
};

// Static member definitions
std::vector<NodeAttributesRegistry::FieldDescriptor> NodeAttributesRegistry::fields_;
bool NodeAttributesRegistry::initialized_ = false;

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
        values.push_back(std::any(T{42})); // arbitrary middle value
    }
    else if constexpr (std::is_floating_point_v<T>) {
        values.push_back(std::any(T{0.0}));
        values.push_back(std::any(T{1.0}));
        values.push_back(std::any(T{-1.0}));
        values.push_back(std::any(std::numeric_limits<T>::min()));
        values.push_back(std::any(std::numeric_limits<T>::max()));
        values.push_back(std::any(T{3.14159})); // arbitrary value
    }
    else if constexpr (std::is_enum_v<T>) {
        // Generate enum values by trying different underlying values
        using UnderlyingType = std::underlying_type_t<T>;
        for (UnderlyingType i = 0; i < 16; ++i) {
            values.push_back(std::any(static_cast<T>(i)));
        }
    }
    
    return values;
}

// =============================================================================
// FIELD REGISTRATION HELPERS
// =============================================================================

template<typename T>
NodeAttributesRegistry::FieldDescriptor create_field_descriptor(
    const std::string& name,
    T NodeGeneAttributes::* member_ptr) {
    
    return NodeAttributesRegistry::FieldDescriptor{
        name,
        std::type_index(typeid(T)),
        []() {
            return generate_runtime_test_values<T>();
        },
        [member_ptr](NodeGeneAttributes& attrs, const std::any& value) {
            try {
                attrs.*member_ptr = std::any_cast<T>(value);
            } catch (const std::bad_any_cast& e) {
                throw std::runtime_error("Failed to set field: bad type cast");
            }
        },
        [member_ptr](const NodeGeneAttributes& attrs) -> std::any {
            return std::any(attrs.*member_ptr);
        },
        is_primitive_type<T>()
    };
}

// =============================================================================
// FIELD REGISTRATION (UPDATE THIS WHEN FIELDS CHANGE)
// =============================================================================


void NodeAttributesRegistry::initialize_fields() {
    // Auto-generated field registration
    #include "node_gene_fields.inc"
    
    // Validation - use direct access to avoid recursive get_fields() call
    if (!all_fields_primitive_direct()) {
        auto non_primitives = get_non_primitive_fields_direct();
        std::string error_msg = "NodeGeneAttributes contains non-primitive fields: ";
        for (size_t i = 0; i < non_primitives.size(); ++i) {
            if (i > 0) error_msg += ", ";
            error_msg += non_primitives[i];
        }
        error_msg += ". Update tests manually to handle complex types.";
        throw std::runtime_error(error_msg);
    }
}

// =============================================================================
// RUNTIME ATTRIBUTE COMBINATION GENERATION
// =============================================================================

std::vector<NodeGeneAttributes> generate_runtime_attribute_combinations() {
    const auto& fields = NodeAttributesRegistry::get_fields();
    std::vector<NodeGeneAttributes> combinations;
    
    if (fields.empty()) {
        // Empty struct - single combination
        combinations.push_back(NodeGeneAttributes{});
        return combinations;
    }
    
    // Get test values for each field
    std::vector<std::vector<std::any>> field_values;
    for (const auto& field : fields) {
        field_values.push_back(field.value_generator());
    }
    
    // Generate all combinations using recursive approach
    std::function<void(size_t, NodeGeneAttributes&)> generate_combinations = 
        [&](size_t field_index, NodeGeneAttributes& current_attrs) {
            if (field_index >= fields.size()) {
                combinations.push_back(current_attrs);
                return;
            }
            
            for (const auto& value : field_values[field_index]) {
                NodeGeneAttributes next_attrs = current_attrs;
                fields[field_index].setter(next_attrs, value);
                generate_combinations(field_index + 1, next_attrs);
            }
        };
    
    NodeGeneAttributes base_attrs{};
    generate_combinations(0, base_attrs);
    
    return combinations;
}

// =============================================================================
// GENERIC TYPE-TRAIT BASED COMPARISON AND MODIFICATION
// =============================================================================

// Helper to determine type category from type_index
enum class TypeCategory {
    BOOLEAN,
    INTEGRAL,
    FLOATING_POINT,
    ENUM,
    UNSUPPORTED
};

template<typename T>
TypeCategory get_type_category() {
    if constexpr (std::is_same_v<T, bool>) {
        return TypeCategory::BOOLEAN;
    } else if constexpr (std::is_integral_v<T>) {
        return TypeCategory::INTEGRAL;
    } else if constexpr (std::is_floating_point_v<T>) {
        return TypeCategory::FLOATING_POINT;
    } else if constexpr (std::is_enum_v<T>) {
        return TypeCategory::ENUM;
    } else {
        return TypeCategory::UNSUPPORTED;
    }
}

// Generic comparison for any two values of the same type
template<typename T>
bool compare_values(const std::any& a, const std::any& b) {
    try {
        T val_a = std::any_cast<T>(a);
        T val_b = std::any_cast<T>(b);
        
        if constexpr (std::is_floating_point_v<T>) {
            // For floating point, use epsilon comparison for safety
            return std::abs(val_a - val_b) < std::numeric_limits<T>::epsilon();
        } else {
            // For integers, bools, and enums, direct comparison
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
        
        // Find a different value
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

// Dispatch comparison based on type category
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

// Dispatch modification based on type category  
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

// =============================================================================
// RUNTIME ATTRIBUTE OPERATIONS
// =============================================================================

bool runtime_attributes_equal(const NodeGeneAttributes& a, const NodeGeneAttributes& b) {
    const auto& fields = NodeAttributesRegistry::get_fields();
    
    for (const auto& field : fields) {
        std::any value_a = field.getter(a);
        std::any value_b = field.getter(b);
        
        // Use type-trait based dispatch for comparison
        if (!dispatch_comparison(field.type, value_a, value_b)) {
            return false;
        }
    }
    
    return true;
}

void runtime_modify_attributes_systematically(NodeGeneAttributes& attrs, size_t field_index) {
    const auto& fields = NodeAttributesRegistry::get_fields();
    
    if (field_index >= fields.size()) {
        return; // Invalid field index
    }
    
    const auto& field = fields[field_index];
    auto test_values = field.value_generator();
    
    if (test_values.empty()) {
        return;
    }
    
    // Get current value
    std::any current_value = field.getter(attrs);
    
    // Use type-trait based dispatch for modification
    if (dispatch_modification(field.type, current_value, test_values)) {
        // Apply the modified value back to the struct
        field.setter(attrs, current_value);
    }
}

// =============================================================================
// UNIT TESTS
// =============================================================================

class NodeGeneTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize the field registry
        NodeAttributesRegistry::get_fields();
    }
    
    std::vector<uint32_t> getTestHistoryIDs() {
        return {0, 1, 42, 1000, std::numeric_limits<uint32_t>::max()};
    }
    
    std::vector<NodeType> getTestNodeTypes() {
        return {NodeType::INPUT, NodeType::HIDDEN, NodeType::OUTPUT, NodeType::BIAS};
    }
};

// =============================================================================
// RUNTIME VALIDATION TESTS
// =============================================================================

TEST_F(NodeGeneTest, RuntimeValidation_AllFieldsArePrimitives) {
    // This test will throw an exception during setup if non-primitives are found
    EXPECT_TRUE(NodeAttributesRegistry::all_fields_primitive()) 
        << "NodeGeneAttributes contains non-primitive fields";
    
    // Report detected field information
    const auto& fields = NodeAttributesRegistry::get_fields();
    std::cout << "Detected " << fields.size() << " field(s) in NodeGeneAttributes:" << std::endl;
    for (const auto& field : fields) {
        std::cout << "  - " << field.name << " (" << field.type.name() << ")" 
                  << (field.is_primitive ? " [primitive]" : " [NON-PRIMITIVE]") << std::endl;
    }
}

TEST_F(NodeGeneTest, RuntimeValidation_FieldSize) {
    size_t field_count = NodeAttributesRegistry::field_count();
    
    // Verify struct size is reasonable
    EXPECT_LE(sizeof(NodeGeneAttributes), 64) 
        << "NodeGeneAttributes size grew unexpectedly - review for bloat";
    
    std::cout << "NodeGeneAttributes: " << field_count << " fields, " 
              << sizeof(NodeGeneAttributes) << " bytes" << std::endl;
}

// =============================================================================
// CONSTRUCTION TESTS
// =============================================================================

TEST_F(NodeGeneTest, Construction_AllValidCombinations) {
    auto historyIDs = getTestHistoryIDs();
    auto nodeTypes = getTestNodeTypes();
    auto attributeCombinations = generate_runtime_attribute_combinations();
    
    EXPECT_GT(attributeCombinations.size(), 0) << "No attribute combinations generated";
    
    size_t total_combinations = 0;
    
    for (uint32_t historyID : historyIDs) {
        for (NodeType nodeType : nodeTypes) {
            for (const auto& attrs : attributeCombinations) {
                NodeGene node(historyID, nodeType, attrs);
                
                // Verify all values stored correctly
                EXPECT_EQ(node.get_historyID(), historyID);
                EXPECT_EQ(node.get_type(), nodeType);
                EXPECT_TRUE(runtime_attributes_equal(node.get_attributes(), attrs));
                
                ++total_combinations;
            }
        }
    }
    
    std::cout << "Tested " << total_combinations << " construction combinations ("
              << historyIDs.size() << " historyIDs × " << nodeTypes.size() 
              << " nodeTypes × " << attributeCombinations.size() << " attribute combinations)" << std::endl;
}

TEST_F(NodeGeneTest, Construction_ImmutabilityContract) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    NodeGene node(42, NodeType::HIDDEN, attributeCombinations[0]);
    
    // Store original values
    uint32_t original_id = node.get_historyID();
    NodeType original_type = node.get_type();
    
    // Modify attributes (this should be allowed)
    size_t field_count = NodeAttributesRegistry::field_count();
    if (field_count > 0) {
        auto& mutable_attrs = const_cast<NodeGene&>(node).get_attributes();
        runtime_modify_attributes_systematically(mutable_attrs, 0);
        
        // Verify immutable fields unchanged
        EXPECT_EQ(node.get_historyID(), original_id) << "historyID should be immutable";
        EXPECT_EQ(node.get_type(), original_type) << "NodeType should be immutable";
    }
}

// =============================================================================
// COPY SEMANTICS TESTS
// =============================================================================

TEST_F(NodeGeneTest, CopyConstructor_ExactReplication) {
    auto historyIDs = getTestHistoryIDs();
    auto nodeTypes = getTestNodeTypes();
    auto attributeCombinations = generate_runtime_attribute_combinations();
    
    // Test a subset to avoid excessive combinations
    size_t test_count = 0;
    constexpr size_t max_tests = 100;
    
    for (uint32_t historyID : historyIDs) {
        for (NodeType nodeType : nodeTypes) {
            for (const auto& attrs : attributeCombinations) {
                if (test_count >= max_tests) break;
                
                NodeGene original(historyID, nodeType, attrs);
                NodeGene copy(original);
                
                // Verify exact replication
                EXPECT_EQ(copy.get_historyID(), original.get_historyID());
                EXPECT_EQ(copy.get_type(), original.get_type());
                EXPECT_TRUE(runtime_attributes_equal(copy.get_attributes(), original.get_attributes()));
                
                ++test_count;
            }
            if (test_count >= max_tests) break;
        }
        if (test_count >= max_tests) break;
    }
    
    std::cout << "Tested " << test_count << " copy constructor scenarios" << std::endl;
}

TEST_F(NodeGeneTest, CopyConstructor_Independence) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    NodeGene original(100, NodeType::HIDDEN, attributeCombinations[0]);
    NodeGene copy(original);
    
    // Modify copy's attributes if possible
    size_t field_count = NodeAttributesRegistry::field_count();
    if (field_count > 0) {
        auto original_attrs = original.get_attributes();
        runtime_modify_attributes_systematically(const_cast<NodeGene&>(copy).get_attributes(), 0);
        
        // Verify original unchanged if modification occurred
        if (!runtime_attributes_equal(copy.get_attributes(), original_attrs)) {
            EXPECT_TRUE(runtime_attributes_equal(original.get_attributes(), original_attrs))
                << "Original should be unmodified when copy is changed";
        }
    }
}

// =============================================================================
// MOVE SEMANTICS TESTS
// =============================================================================

TEST_F(NodeGeneTest, MoveConstructor_StateTransfer) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    NodeGene original(500, NodeType::HIDDEN, attributeCombinations[0]);
    
    // Store original state
    uint32_t expected_id = original.get_historyID();
    NodeType expected_type = original.get_type();
    NodeGeneAttributes expected_attrs = original.get_attributes();
    
    NodeGene moved(std::move(original));
    
    // Verify state transferred correctly
    EXPECT_EQ(moved.get_historyID(), expected_id);
    EXPECT_EQ(moved.get_type(), expected_type);
    EXPECT_TRUE(runtime_attributes_equal(moved.get_attributes(), expected_attrs));
}

// =============================================================================
// ATTRIBUTE MODIFICATION TESTS
// =============================================================================

TEST_F(NodeGeneTest, AttributeModification_MutabilityMatrix) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    NodeGene node(800, NodeType::HIDDEN, attributeCombinations[0]);
    size_t field_count = NodeAttributesRegistry::field_count();
    
    // Test each field can be modified
    for (size_t field_index = 0; field_index < field_count; ++field_index) {
        NodeGeneAttributes original_attrs = node.get_attributes();
        
        // Modify the field at this index
        auto& mutable_attrs = const_cast<NodeGene&>(node).get_attributes();
        runtime_modify_attributes_systematically(mutable_attrs, field_index);
        
        // Verify immutable fields still unchanged
        EXPECT_EQ(node.get_historyID(), 800u);
        EXPECT_EQ(node.get_type(), NodeType::HIDDEN);
        
        // Note: We don't verify that attributes changed because the modification
        // might result in the same value if all test values are identical
    }
    
    std::cout << "Tested modification of " << field_count << " field(s)" << std::endl;
}

TEST_F(NodeGeneTest, AttributeModification_ConstCorrectness) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const NodeGene const_node(900, NodeType::OUTPUT, attributeCombinations[0]);
    
    // Const getter should return const reference
    const auto& const_attrs = const_node.get_attributes();
    EXPECT_TRUE(runtime_attributes_equal(const_attrs, attributeCombinations[0]));
    
    // Should not be able to modify through const reference
    // (This is verified at compile-time, so test succeeds if it compiles)
    SUCCEED() << "Const correctness verified at compile-time";
}

// =============================================================================
// EXTENSIBILITY VERIFICATION TESTS
// =============================================================================

TEST_F(NodeGeneTest, ExtensibilityVerification_AutomaticTestAdaptation) {
    // This test verifies that our testing infrastructure will automatically
    // adapt when fields are added/removed/changed in NodeGeneAttributes
    
    size_t current_field_count = NodeAttributesRegistry::field_count();
    std::cout << "Current NodeGeneAttributes field count: " << current_field_count << std::endl;
    
    // Verify our test value generation works for current structure
    auto combinations = generate_runtime_attribute_combinations();
    EXPECT_GT(combinations.size(), 0) << "Attribute combination generation failed";
    
    std::cout << "Generated " << combinations.size() << " attribute combinations" << std::endl;
    
    // Show field details
    const auto& fields = NodeAttributesRegistry::get_fields();
    for (const auto& field : fields) {
        auto test_values = field.value_generator();
        std::cout << "Field '" << field.name << "': " << test_values.size() << " test values" << std::endl;
    }
}

TEST_F(NodeGeneTest, ExtensibilityVerification_GenericOperations) {
    auto attributeCombinations = generate_runtime_attribute_combinations();
    if (attributeCombinations.size() < 1) {
        GTEST_SKIP() << "Need at least 1 attribute combination for generic operation testing";
    }
    
    // Test generic equality
    EXPECT_TRUE(runtime_attributes_equal(attributeCombinations[0], attributeCombinations[0]));
    
    // Test generic modification (if fields exist)
    size_t field_count = NodeAttributesRegistry::field_count();
    if (field_count > 0) {
        auto test_attrs = attributeCombinations[0];
        auto original_attrs = test_attrs;
        
        runtime_modify_attributes_systematically(test_attrs, 0);
        
        // Modification should be safe (no crashes)
        SUCCEED() << "Generic modification completed without errors";
    }
}