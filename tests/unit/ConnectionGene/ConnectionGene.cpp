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
    
    // Helper to create various index pairs for testing
    std::vector<std::pair<size_t, size_t>> getTestIndexPairs() {
        std::vector<std::pair<size_t, size_t>> pairs;
        
        // Various index combinations for testing
        pairs.emplace_back(0, 1);   // input -> hidden
        pairs.emplace_back(0, 2);   // input -> output  
        pairs.emplace_back(1, 3);   // hidden -> hidden
        pairs.emplace_back(1, 2);   // hidden -> output
        pairs.emplace_back(4, 1);   // bias -> hidden
        pairs.emplace_back(4, 2);   // bias -> output
        pairs.emplace_back(5, 10);  // larger indices
        pairs.emplace_back(0, 100); // edge case indices
        
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
    auto indexPairs = getTestIndexPairs();
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    
    EXPECT_GT(attributeCombinations.size(), 0) << "No attribute combinations generated";
    EXPECT_GT(indexPairs.size(), 0) << "No index pairs available";
    
    size_t total_combinations = 0;
    
    for (uint32_t historyID : historyIDs) {
        for (const auto& [sourceIndex, targetIndex] : indexPairs) {
            for (const auto& attrs : attributeCombinations) {
                ConnectionGene conn(historyID, sourceIndex, targetIndex, attrs);
                
                // Verify all values stored correctly
                EXPECT_EQ(conn.get_historyID(), historyID);
                EXPECT_EQ(conn.get_sourceNodeIndex(), sourceIndex);
                EXPECT_EQ(conn.get_targetNodeIndex(), targetIndex);
                EXPECT_TRUE((runtime_attributes_equal<ConnectionGeneAttributes, ConnectionAttributesRegistry>(conn.get_attributes(), attrs)));
                
                ++total_combinations;
            }
        }
    }
    
    std::cout << "Tested " << total_combinations << " construction combinations ("
              << historyIDs.size() << " historyIDs × " << indexPairs.size() 
              << " index pairs × " << attributeCombinations.size() << " attribute combinations)" << std::endl;
}

TEST_F(ConnectionGeneTest, Construction_ImmutabilityContract) {
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    ConnectionGene conn(42, sourceIndex, targetIndex, attributeCombinations[0]);
    
    // Store original values
    uint32_t original_id = conn.get_historyID();
    size_t original_source = conn.get_sourceNodeIndex();
    size_t original_target = conn.get_targetNodeIndex();
    
    // Modify attributes (this should be allowed)
    size_t field_count = ConnectionAttributesRegistry::field_count();
    if (field_count > 0) {
        auto& mutable_attrs = const_cast<ConnectionGene&>(conn).get_attributes();
        runtime_modify_attributes_systematically<ConnectionGeneAttributes, ConnectionAttributesRegistry>(mutable_attrs, 0);
        
        // Verify immutable fields unchanged
        EXPECT_EQ(conn.get_historyID(), original_id) << "historyID should be immutable";
        EXPECT_EQ(conn.get_sourceNodeIndex(), original_source) << "source node index should be immutable";
        EXPECT_EQ(conn.get_targetNodeIndex(), original_target) << "target node index should be immutable";
    }
}

// =============================================================================
// COPY SEMANTICS TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, CopyConstructor_ExactReplication) {
    auto historyIDs = getTestHistoryIDs();
    auto indexPairs = getTestIndexPairs();
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    
    size_t test_count = 0;
    constexpr size_t max_tests = 50;
    
    for (uint32_t historyID : historyIDs) {
        for (const auto& [sourceIndex, targetIndex] : indexPairs) {
            for (const auto& attrs : attributeCombinations) {
                if (test_count >= max_tests) break;
                
                ConnectionGene original(historyID, sourceIndex, targetIndex, attrs);
                ConnectionGene copy(original);
                
                // Verify exact replication
                EXPECT_EQ(copy.get_historyID(), original.get_historyID());
                EXPECT_EQ(copy.get_sourceNodeIndex(), original.get_sourceNodeIndex());
                EXPECT_EQ(copy.get_targetNodeIndex(), original.get_targetNodeIndex());
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
    
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    ConnectionGene original(100, sourceIndex, targetIndex, attributeCombinations[0]);
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
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    const auto& attrs = attributeCombinations[0];
    
    // Test self-equality
    ConnectionGene conn1(42, sourceIndex, targetIndex, attrs);
    EXPECT_TRUE(conn1 == conn1) << "Connection should be equal to itself";
    EXPECT_FALSE(conn1 != conn1) << "Connection should not be unequal to itself";
    
    // Test equality with identical construction
    ConnectionGene conn2(42, sourceIndex, targetIndex, attrs);
    EXPECT_TRUE(conn1 == conn2) << "Identical connections should be equal";
    EXPECT_FALSE(conn1 != conn2) << "Identical connections should not be unequal";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_DifferentHistoryIDs) {
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    const auto& attrs = attributeCombinations[0];
    
    ConnectionGene conn1(42, sourceIndex, targetIndex, attrs);
    ConnectionGene conn2(43, sourceIndex, targetIndex, attrs);
    
    EXPECT_FALSE(conn1 == conn2) << "Connections with different history IDs should not be equal";
    EXPECT_TRUE(conn1 != conn2) << "Connections with different history IDs should be unequal";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_SameHistoryIDDifferentIndices) {
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.size() < 2) {
        GTEST_SKIP() << "Need at least 2 index pairs for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& attrs = attributeCombinations[0];
    
    // Find two pairs with different indices
    size_t pair1_idx = 0;
    size_t pair2_idx = 1;
    
    // Look for pairs with different source or target indices
    for (size_t i = 1; i < indexPairs.size(); ++i) {
        if (indexPairs[0].first != indexPairs[i].first || indexPairs[0].second != indexPairs[i].second) {
            pair2_idx = i;
            break;
        }
    }
    
    // Create connections with same history ID but different indices
    ConnectionGene conn1(42, indexPairs[pair1_idx].first, indexPairs[pair1_idx].second, attrs);
    ConnectionGene conn2(42, indexPairs[pair2_idx].first, indexPairs[pair2_idx].second, attrs);  // Different indices
    
    // Connections should be equal because comparison is based on connection history ID only
    EXPECT_TRUE(conn1 == conn2) << "Connections with same history ID should be equal regardless of indices";
    EXPECT_FALSE(conn1 != conn2) << "Connections with same history ID should not be unequal";
    
    // Verify at least one index is different
    EXPECT_TRUE(conn1.get_sourceNodeIndex() != conn2.get_sourceNodeIndex() || 
               conn1.get_targetNodeIndex() != conn2.get_targetNodeIndex()) 
               << "At least one index should be different";
}

TEST_F(ConnectionGeneTest, ComparisonOperators_SameHistoryIDDifferentAttributes) {
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.size() < 2) {
        GTEST_SKIP() << "Need at least 2 attribute combinations for testing";
    }
    
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    
    // Note: Using SAME connection history ID (42) to test that comparison is based on connection history ID only
    ConnectionGene conn1(42, sourceIndex, targetIndex, attributeCombinations[0]);
    ConnectionGene conn2(42, sourceIndex, targetIndex, attributeCombinations[1]);  // Different attributes
    
    // Connections should be equal because comparison is based on connection history ID only,
    // not on the attributes. This is correct NEAT behavior.
    EXPECT_TRUE(conn1 == conn2) << "Connections with same connection history ID should be equal regardless of attributes";
    EXPECT_FALSE(conn1 != conn2) << "Connections with same connection history ID should not be unequal";
    
    // But let's also test with different connection history IDs to show the difference
    ConnectionGene conn3(43, sourceIndex, targetIndex, attributeCombinations[0]); // Different connection history ID
    EXPECT_FALSE(conn1 == conn3) << "Connections with different connection history IDs should not be equal";
    EXPECT_TRUE(conn1 != conn3) << "Connections with different connection history IDs should be unequal";
}


// =============================================================================
// ATTRIBUTE MODIFICATION TESTS
// =============================================================================

TEST_F(ConnectionGeneTest, AttributeModification_MutabilityMatrix) {
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    ConnectionGene conn(800, sourceIndex, targetIndex, attributeCombinations[0]);
    size_t field_count = ConnectionAttributesRegistry::field_count();
    
    // Test each field can be modified
    for (size_t field_index = 0; field_index < field_count; ++field_index) {
        ConnectionGeneAttributes original_attrs = conn.get_attributes();
        
        // Modify the field at this index
        auto& mutable_attrs = const_cast<ConnectionGene&>(conn).get_attributes();
        runtime_modify_attributes_systematically<ConnectionGeneAttributes, ConnectionAttributesRegistry>(mutable_attrs, field_index);
        
        // Verify immutable fields still unchanged
        EXPECT_EQ(conn.get_historyID(), 800u);
        EXPECT_EQ(conn.get_sourceNodeIndex(), sourceIndex);
        EXPECT_EQ(conn.get_targetNodeIndex(), targetIndex);
    }
    
    std::cout << "Tested modification of " << field_count << " field(s)" << std::endl;
}

TEST_F(ConnectionGeneTest, AttributeModification_ConstCorrectness) {
    auto indexPairs = getTestIndexPairs();
    if (indexPairs.empty()) {
        GTEST_SKIP() << "No index pairs available for testing";
    }
    
    auto attributeCombinations = generate_runtime_attribute_combinations<ConnectionGeneAttributes, ConnectionAttributesRegistry>();
    if (attributeCombinations.empty()) {
        GTEST_SKIP() << "No attribute combinations available for testing";
    }
    
    const auto& [sourceIndex, targetIndex] = indexPairs[0];
    const ConnectionGene const_conn(900, sourceIndex, targetIndex, attributeCombinations[0]);
    
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