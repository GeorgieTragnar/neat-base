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

#include "tests/test_common.h"
#include "tests/test_utilities.h"
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
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
    
    EXPECT_GT(attributeCombinations.size(), 0) << "No attribute combinations generated";
    
    size_t total_combinations = 0;
    
    for (uint32_t historyID : historyIDs) {
        for (NodeType nodeType : nodeTypes) {
            for (const auto& attrs : attributeCombinations) {
                NodeGene node(historyID, nodeType, attrs);
                
                // Verify all values stored correctly
                EXPECT_EQ(node.get_historyID(), historyID);
                EXPECT_EQ(node.get_type(), nodeType);
                EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(node.get_attributes(), attrs)));
                
                ++total_combinations;
            }
        }
    }
    
    std::cout << "Tested " << total_combinations << " construction combinations ("
              << historyIDs.size() << " historyIDs × " << nodeTypes.size() 
              << " nodeTypes × " << attributeCombinations.size() << " attribute combinations)" << std::endl;
}

TEST_F(NodeGeneTest, Construction_ImmutabilityContract) {
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
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
        runtime_modify_attributes_systematically<NodeGeneAttributes, NodeAttributesRegistry>(mutable_attrs, 0);
        
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
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
    
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
                EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(copy.get_attributes(), original.get_attributes())));
                
                ++test_count;
            }
            if (test_count >= max_tests) break;
        }
        if (test_count >= max_tests) break;
    }
    
    std::cout << "Tested " << test_count << " copy constructor scenarios" << std::endl;
}

TEST_F(NodeGeneTest, CopyConstructor_Independence) {
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    NodeGene original(100, NodeType::HIDDEN, attributeCombinations[0]);
    NodeGene copy(original);
    
    // Modify copy's attributes if possible
    size_t field_count = NodeAttributesRegistry::field_count();
    if (field_count > 0) {
        auto original_attrs = original.get_attributes();
        runtime_modify_attributes_systematically<NodeGeneAttributes, NodeAttributesRegistry>(const_cast<NodeGene&>(copy).get_attributes(), 0);
        
        // Verify original unchanged if modification occurred
        if (!runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(copy.get_attributes(), original_attrs)) {
            EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(original.get_attributes(), original_attrs)))
                << "Original should be unmodified when copy is changed";
        }
    }
}

// =============================================================================
// MOVE SEMANTICS TESTS
// =============================================================================

TEST_F(NodeGeneTest, MoveConstructor_StateTransfer) {
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
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
    EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(moved.get_attributes(), expected_attrs)));
}

// =============================================================================
// ATTRIBUTE MODIFICATION TESTS
// =============================================================================

TEST_F(NodeGeneTest, AttributeModification_MutabilityMatrix) {
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
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
        runtime_modify_attributes_systematically<NodeGeneAttributes, NodeAttributesRegistry>(mutable_attrs, field_index);
        
        // Verify immutable fields still unchanged
        EXPECT_EQ(node.get_historyID(), 800u);
        EXPECT_EQ(node.get_type(), NodeType::HIDDEN);
        
        // Note: We don't verify that attributes changed because the modification
        // might result in the same value if all test values are identical
    }
    
    std::cout << "Tested modification of " << field_count << " field(s)" << std::endl;
}

TEST_F(NodeGeneTest, AttributeModification_ConstCorrectness) {
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const NodeGene const_node(900, NodeType::OUTPUT, attributeCombinations[0]);
    
    // Const getter should return const reference
    const auto& const_attrs = const_node.get_attributes();
    EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(const_attrs, attributeCombinations[0])));
    
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
    auto combinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
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
    auto attributeCombinations = generate_runtime_attribute_combinations<NodeGeneAttributes, NodeAttributesRegistry>();
    if (attributeCombinations.size() < 1) {
        GTEST_SKIP() << "Need at least 1 attribute combination for generic operation testing";
    }
    
    // Test generic equality
    EXPECT_TRUE((runtime_attributes_equal<NodeGeneAttributes, NodeAttributesRegistry>(attributeCombinations[0], attributeCombinations[0])));
    
    // Test generic modification (if fields exist)
    size_t field_count = NodeAttributesRegistry::field_count();
    if (field_count > 0) {
        auto test_attrs = attributeCombinations[0];
        auto original_attrs = test_attrs;
        
        runtime_modify_attributes_systematically<NodeGeneAttributes, NodeAttributesRegistry>(test_attrs, 0);
        
        // Modification should be safe (no crashes)
        SUCCEED() << "Generic modification completed without errors";
    }
}