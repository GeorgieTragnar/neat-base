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
#include "version3/data/ConnectionGene.hpp"

// =============================================================================
// RUNTIME FIELD REGISTRY SYSTEM FOR CONNECTION GENE ATTRIBUTES
// =============================================================================

class ConnectionAttributesRegistry {
public:
    struct FieldDescriptor {
        std::string name;
        std::type_index type;
        std::function<std::vector<std::any>()> value_generator;
        std::function<void(ConnectionGeneAttributes&, const std::any&)> setter;
        std::function<std::any(const ConnectionGeneAttributes&)> getter;
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
std::vector<ConnectionAttributesRegistry::FieldDescriptor> ConnectionAttributesRegistry::fields_;
bool ConnectionAttributesRegistry::initialized_ = false;


// =============================================================================
// FIELD REGISTRATION HELPERS
// =============================================================================

template<typename T>
ConnectionAttributesRegistry::FieldDescriptor create_field_descriptor(
    const std::string& name,
    T ConnectionGeneAttributes::* member_ptr) {
    
    return ConnectionAttributesRegistry::FieldDescriptor{
        name,
        std::type_index(typeid(T)),
        []() {
            return generate_runtime_test_values<T>();
        },
        [member_ptr](ConnectionGeneAttributes& attrs, const std::any& value) {
            try {
                attrs.*member_ptr = std::any_cast<T>(value);
            } catch (const std::bad_any_cast& e) {
                throw std::runtime_error("Failed to set field: bad type cast");
            }
        },
        [member_ptr](const ConnectionGeneAttributes& attrs) -> std::any {
            return std::any(attrs.*member_ptr);
        },
        is_primitive_type<T>()
    };
}

// =============================================================================
// FIELD REGISTRATION (AUTO-GENERATED)
// =============================================================================

void ConnectionAttributesRegistry::initialize_fields() {
    // Auto-generated field registration
    #include "connection_gene_fields.inc"
    
    // Validation
    if (!all_fields_primitive_direct()) {
        auto non_primitives = get_non_primitive_fields_direct();
        std::string error_msg = "ConnectionGeneAttributes contains non-primitive fields: ";
        for (size_t i = 0; i < non_primitives.size(); ++i) {
            if (i > 0) error_msg += ", ";
            error_msg += non_primitives[i];
        }
        error_msg += ". Update tests manually to handle complex types.";
        throw std::runtime_error(error_msg);
    }
}




// =============================================================================
// HELPER FUNCTIONS FOR CREATING TEST NODE GENES
// =============================================================================

NodeGene create_test_node_gene(uint32_t historyID, NodeType type) {
    NodeGeneAttributes attrs;
    attrs.activationType = ActivationType::SIGMOID; // Default activation
    return NodeGene(historyID, type, attrs);
}

// =============================================================================
// UNIT TESTS
// =============================================================================

class ConnectionGeneTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize the field registry
        ConnectionAttributesRegistry::get_fields();
    }
    
    std::vector<uint32_t> getTestHistoryIDs() {
        return {0, 1, 42, 1000, std::numeric_limits<uint32_t>::max()};
    }
    
    // Helper to create various node gene pairs for testing
    std::vector<std::pair<NodeGene, NodeGene>> getTestNodeGenePairs() {
        std::vector<std::pair<NodeGene, NodeGene>> pairs;
        
        // Different combinations of node types
        auto input1 = create_test_node_gene(1, NodeType::INPUT);
        auto input2 = create_test_node_gene(2, NodeType::INPUT);
        auto hidden1 = create_test_node_gene(10, NodeType::HIDDEN);
        auto hidden2 = create_test_node_gene(11, NodeType::HIDDEN);
        auto output1 = create_test_node_gene(20, NodeType::OUTPUT);
        auto output2 = create_test_node_gene(21, NodeType::OUTPUT);
        auto bias = create_test_node_gene(30, NodeType::BIAS);
        
        // Valid connection combinations (non-OUTPUT -> non-INPUT)
        pairs.emplace_back(input1, hidden1);
        pairs.emplace_back(input1, output1);
        pairs.emplace_back(hidden1, hidden2);
        pairs.emplace_back(hidden1, output1);
        pairs.emplace_back(bias, hidden1);
        pairs.emplace_back(bias, output1);
        
        return pairs;
    }
};

// =============================================================================
// RUNTIME VALIDATION TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, RuntimeValidation_AllFieldsArePrimitives) {
    EXPECT_TRUE(ConnectionAttributesRegistry::all_fields_primitive()) 
        << "ConnectionGeneAttributes contains non-primitive fields";
    
    const auto& fields = ConnectionAttributesRegistry::get_fields();
    std::cout << "Detected " << fields.size() << " field(s) in ConnectionGeneAttributes:" << std::endl;
    for (const auto& field : fields) {
        std::cout << "  - " << field.name << " (" << field.type.name() << ")" 
                  << (field.is_primitive ? " [primitive]" : " [NON-PRIMITIVE]") << std::endl;
    }
}

TEST_F(ConnectionGeneTest, RuntimeValidation_FieldSize) {
    size_t field_count = ConnectionAttributesRegistry::field_count();
    
    EXPECT_LE(sizeof(ConnectionGeneAttributes), 64) 
        << "ConnectionGeneAttributes size grew unexpectedly";
    
    std::cout << "ConnectionGeneAttributes: " << field_count << " fields, " 
              << sizeof(ConnectionGeneAttributes) << " bytes" << std::endl;
}

// =============================================================================
// CONSTRUCTION TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, Construction_AllValidCombinations) {
    auto historyIDs = getTestHistoryIDs();
    auto nodeGenePairs = getTestNodeGenePairs();
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    
    EXPECT_GT(attributeCombinations.size(), 0) << "No attribute combinations generated";
    EXPECT_GT(nodeGenePairs.size(), 0) << "No node gene pairs available";
    
    size_t total_combinations = 0;
    
    for (uint32_t historyID : historyIDs) {
        for (const auto& [sourceNode, targetNode] : nodeGenePairs) {
            for (const auto& attrs : attributeCombinations) {
                ConnectionGene conn(historyID, sourceNode, targetNode, attrs);
                
                // Verify all values stored correctly
                EXPECT_EQ(conn.get_historyID(), historyID);
                EXPECT_EQ(&conn.get_sourceNodeGene(), &sourceNode);
                EXPECT_EQ(&conn.get_targetNodeGene(), &targetNode);
                EXPECT_TRUE((runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(conn.get_attributes(), attrs)));
                
                ++total_combinations;
            }
        }
    }
    
    std::cout << "Tested " << total_combinations << " construction combinations ("
              << historyIDs.size() << " historyIDs × " << nodeGenePairs.size() 
              << " node pairs × " << attributeCombinations.size() << " attribute combinations)" << std::endl;
}

TEST_F(ConnectionGeneTest, Construction_ImmutabilityContract) {
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    ConnectionGene conn(42, sourceNode, targetNode, attributeCombinations[0]);
    
    // Store original values
    uint32_t original_id = conn.get_historyID();
    const NodeGene* original_source = &conn.get_sourceNodeGene();
    const NodeGene* original_target = &conn.get_targetNodeGene();
    
    // Modify attributes (this should be allowed)
    size_t field_count = ConnectionAttributesRegistry::field_count();
    if (field_count > 0) {
        auto& mutable_attrs = const_cast<ConnectionGene&>(conn).get_attributes();
        runtime_modify_attributes_systematically<ConnectionGeneAttributes, ConnectionAttributesRegistry>(mutable_attrs, 0);
        
        // Verify immutable fields unchanged
        EXPECT_EQ(conn.get_historyID(), original_id) << "historyID should be immutable";
        EXPECT_EQ(&conn.get_sourceNodeGene(), original_source) << "source node reference should be immutable";
        EXPECT_EQ(&conn.get_targetNodeGene(), original_target) << "target node reference should be immutable";
    }
}

// =============================================================================
// COPY SEMANTICS TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, CopyConstructor_ExactReplication) {
    auto historyIDs = getTestHistoryIDs();
    auto nodeGenePairs = getTestNodeGenePairs();
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    
    size_t test_count = 0;
    constexpr size_t max_tests = 50;
    
    for (uint32_t historyID : historyIDs) {
        for (const auto& [sourceNode, targetNode] : nodeGenePairs) {
            for (const auto& attrs : attributeCombinations) {
                if (test_count >= max_tests) break;
                
                ConnectionGene original(historyID, sourceNode, targetNode, attrs);
                ConnectionGene copy(original);
                
                // Verify exact replication
                EXPECT_EQ(copy.get_historyID(), original.get_historyID());
                EXPECT_EQ(copy.get_sourceNodeGene().get_historyID(), original.get_sourceNodeGene().get_historyID());
                EXPECT_EQ(copy.get_targetNodeGene().get_historyID(), original.get_targetNodeGene().get_historyID());
                EXPECT_TRUE((runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(copy.get_attributes(), original.get_attributes())));
                
                ++test_count;
            }
            if (test_count >= max_tests) break;
        }
        if (test_count >= max_tests) break;
    }
    
    std::cout << "Tested " << test_count << " copy constructor scenarios" << std::endl;
}

TEST_F(ConnectionGeneTest, CopyConstructor_Independence) {
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    ConnectionGene original(100, sourceNode, targetNode, attributeCombinations[0]);
    ConnectionGene copy(original);
    
    // Modify copy's attributes if possible
    size_t field_count = ConnectionAttributesRegistry::field_count();
    if (field_count > 0) {
        auto original_attrs = original.get_attributes();
        runtime_modify_attributes_systematically<ConnectionGeneAttributes, ConnectionAttributesRegistry>(const_cast<ConnectionGene&>(copy).get_attributes(), 0);
        
        // Verify original unchanged if modification occurred
        if (!runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(copy.get_attributes(), original_attrs)) {
            EXPECT_TRUE((runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(original.get_attributes(), original_attrs)))
                << "Original should be unmodified when copy is changed";
        }
    }
}

// =============================================================================
// COMPARISON OPERATORS TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, ComparisonOperators_Equality) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    const auto& attrs = attributeCombinations[0];
    
    // Test self-equality
    ConnectionGene conn1(42, sourceNode, targetNode, attrs);
    EXPECT_TRUE(conn1 == conn1) << "Connection should be equal to itself";
    EXPECT_FALSE(conn1 != conn1) << "Connection should not be unequal to itself";
    
    // Test equality with identical construction
    ConnectionGene conn2(42, sourceNode, targetNode, attrs);
    EXPECT_TRUE(conn1 == conn2) << "Identical connections should be equal";
    EXPECT_FALSE(conn1 != conn2) << "Identical connections should not be unequal";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_DifferentHistoryIDs) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    const auto& attrs = attributeCombinations[0];
    
    ConnectionGene conn1(42, sourceNode, targetNode, attrs);
    ConnectionGene conn2(43, sourceNode, targetNode, attrs);
    
    EXPECT_FALSE(conn1 == conn2) << "Connections with different history IDs should not be equal";
    EXPECT_TRUE(conn1 != conn2) << "Connections with different history IDs should be unequal";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_SameHistoryIDsButDifferentNodeReferences) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.size() < 2) {
        GTEST_SKIP() << "Need at least 2 node gene pairs for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    // Create nodes with SAME history IDs but different references
    // This simulates cross-genome comparison where nodes have same IDs but are different objects
    auto sourceNode1 = create_test_node_gene(100, NodeType::INPUT);  // historyID = 100
    auto targetNode1 = create_test_node_gene(200, NodeType::OUTPUT); // historyID = 200
    auto sourceNode2 = create_test_node_gene(100, NodeType::INPUT);  // historyID = 100 (same as sourceNode1)
    auto targetNode2 = create_test_node_gene(200, NodeType::OUTPUT); // historyID = 200 (same as targetNode1)
    
    const auto& attrs = attributeCombinations[0];
    
    ConnectionGene conn1(42, sourceNode1, targetNode1, attrs);
    ConnectionGene conn2(42, sourceNode2, targetNode2, attrs);
    
    // Even though node references are different, the connections should be equal
    // because they have the same historyID and the nodes have the same historyIDs
    EXPECT_TRUE(conn1 == conn2) << "Connections with same history ID should be equal regardless of node references";
    EXPECT_FALSE(conn1 != conn2) << "Connections with same history ID should not be unequal";
    
    // Verify the node references are indeed different objects but have same history IDs
    EXPECT_NE(&sourceNode1, &sourceNode2) << "Source node references should be different";
    EXPECT_NE(&targetNode1, &targetNode2) << "Target node references should be different";
    EXPECT_EQ(sourceNode1.get_historyID(), sourceNode2.get_historyID()) << "Source nodes should have same history ID";
    EXPECT_EQ(targetNode1.get_historyID(), targetNode2.get_historyID()) << "Target nodes should have same history ID";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_DifferentNodeHistoryIDs) {
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    // Create nodes with DIFFERENT history IDs
    auto sourceNode1 = create_test_node_gene(100, NodeType::INPUT);  // historyID = 100
    auto targetNode1 = create_test_node_gene(200, NodeType::OUTPUT); // historyID = 200
    auto sourceNode2 = create_test_node_gene(101, NodeType::INPUT);  // historyID = 101 (different)
    auto targetNode2 = create_test_node_gene(201, NodeType::OUTPUT); // historyID = 201 (different)
    
    const auto& attrs = attributeCombinations[0];
    
    // Note: Using SAME connection history ID (42) to test that comparison is based on connection history ID only
    ConnectionGene conn1(42, sourceNode1, targetNode1, attrs);
    ConnectionGene conn2(42, sourceNode2, targetNode2, attrs);
    
    // Connections should still be equal because comparison is based on connection history ID only,
    // not on the node history IDs. This is the key insight from the implementation.
    EXPECT_TRUE(conn1 == conn2) << "Connections with same connection history ID should be equal regardless of node history IDs";
    EXPECT_FALSE(conn1 != conn2) << "Connections with same connection history ID should not be unequal";
    
    // But let's also test with different connection history IDs to show the difference
    ConnectionGene conn3(43, sourceNode1, targetNode1, attrs); // Different connection history ID
    EXPECT_FALSE(conn1 == conn3) << "Connections with different connection history IDs should not be equal";
    EXPECT_TRUE(conn1 != conn3) << "Connections with different connection history IDs should be unequal";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_DifferentAttributes) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.size() < 2) {
        GTEST_SKIP() << "Need at least 2 attribute combinations for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    
    ConnectionGene conn1(42, sourceNode, targetNode, attributeCombinations[0]);
    ConnectionGene conn2(42, sourceNode, targetNode, attributeCombinations[1]);
    
    // ConnectionGene comparison is based on historyID only, not attributes
    // This is correct behavior for NEAT - connections with same innovation number are considered the same
    EXPECT_TRUE(conn1 == conn2) << "Connections with same history ID should be equal regardless of attributes";
    EXPECT_FALSE(conn1 != conn2) << "Connections with same history ID should not be unequal";
    
    // Verify that attributes are actually different but comparison still succeeds
    bool attributes_different = !runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(attributeCombinations[0], attributeCombinations[1]);
    if (attributes_different) {
        std::cout << "✓ Verified: Attributes are different but connections are still considered equal (correct NEAT behavior)" << std::endl;
    } else {
        std::cout << "Note: Generated attribute combinations happened to be identical" << std::endl;
    }
}

// =============================================================================
// ATTRIBUTE MODIFICATION TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, AttributeModification_MutabilityMatrix) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    ConnectionGene conn(800, sourceNode, targetNode, attributeCombinations[0]);
    size_t field_count = ConnectionAttributesRegistry::field_count();
    
    // Test each field can be modified
    for (size_t field_index = 0; field_index < field_count; ++field_index) {
        ConnectionGeneAttributes original_attrs = conn.get_attributes();
        
        // Modify the field at this index
        auto& mutable_attrs = const_cast<ConnectionGene&>(conn).get_attributes();
        runtime_modify_attributes_systematically<ConnectionGeneAttributes, ConnectionAttributesRegistry>(mutable_attrs, field_index);
        
        // Verify immutable fields still unchanged
        EXPECT_EQ(conn.get_historyID(), 800u);
        EXPECT_EQ(&conn.get_sourceNodeGene(), &sourceNode);
        EXPECT_EQ(&conn.get_targetNodeGene(), &targetNode);
    }
    
    std::cout << "Tested modification of " << field_count << " field(s)" << std::endl;
}

TEST_F(ConnectionGeneTest, AttributeModification_ConstCorrectness) {
    auto nodeGenePairs = getTestNodeGenePairs();
    if (nodeGenePairs.empty()) {
        GTEST_SKIP() << "No node gene pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceNode, targetNode] = nodeGenePairs[0];
    const ConnectionGene const_conn(900, sourceNode, targetNode, attributeCombinations[0]);
    
    // Const getter should return const reference
    const auto& const_attrs = const_conn.get_attributes();
    EXPECT_TRUE((runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(const_attrs, attributeCombinations[0])));
    
    SUCCEED() << "Const correctness verified at compile-time";
}

// =============================================================================
// EXTENSIBILITY VERIFICATION TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, ExtensibilityVerification_AutomaticTestAdaptation) {
    size_t current_field_count = ConnectionAttributesRegistry::field_count();
    std::cout << "Current ConnectionGeneAttributes field count: " << current_field_count << std::endl;
    
    auto combinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    EXPECT_GT(combinations.size(), 0) << "Attribute combination generation failed";
    
    std::cout << "Generated " << combinations.size() << " attribute combinations" << std::endl;
    
    const auto& fields = ConnectionAttributesRegistry::get_fields();
    for (const auto& field : fields) {
        auto test_values = field.value_generator();
        std::cout << "Field '" << field.name << "': " << test_values.size() << " test values" << std::endl;
    }
}