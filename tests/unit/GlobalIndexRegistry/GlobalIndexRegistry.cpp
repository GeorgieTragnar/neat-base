#include <gtest/gtest.h>
#include "version3/population/GlobalIndexRegistry.hpp"

using namespace Population;

class GlobalIndexRegistryTest : public ::testing::Test {
protected:
    void SetUp() override {
        registry = std::make_unique<GlobalIndexRegistry>(5);
    }
    
    std::unique_ptr<GlobalIndexRegistry> registry;
};

// ============================================================================
// State Machine Integrity Tests
// ============================================================================

TEST_F(GlobalIndexRegistryTest, StateTransitions_ValidPaths) {
    // Test: All legal state transitions work correctly
    
    // Active -> Elite
    EXPECT_EQ(registry->getState(0), GenomeState::Active);
    registry->markAsElite(0);
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);
    
    // Active -> HotElimination
    EXPECT_EQ(registry->getState(1), GenomeState::Active);
    registry->markForElimination(1);
    EXPECT_EQ(registry->getState(1), GenomeState::HotElimination);
    
    // Elite -> HotElimination
    registry->markForElimination(0);
    EXPECT_EQ(registry->getState(0), GenomeState::HotElimination);
    
    // HotElimination -> ColdElimination
    registry->transitionToCold(0);
    EXPECT_EQ(registry->getState(0), GenomeState::ColdElimination);
    registry->transitionToCold(1);
    EXPECT_EQ(registry->getState(1), GenomeState::ColdElimination);
    
    // ColdElimination -> ReadyForReplacement
    registry->markReadyForReplacement(0);
    EXPECT_EQ(registry->getState(0), GenomeState::ReadyForReplacement);
    registry->markReadyForReplacement(1);
    EXPECT_EQ(registry->getState(1), GenomeState::ReadyForReplacement);
}

TEST_F(GlobalIndexRegistryTest, StateTransitions_InvalidPaths) {
    // Test: Invalid state transitions should assert/fail
    
    // Cannot mark non-Active/Elite for elimination
    registry->markForElimination(0); // Active -> HotElimination (valid)
    EXPECT_DEATH(registry->markForElimination(0), "Can only mark Active or Elite genomes for elimination");
    
    // Cannot transition to cold from non-HotElimination
    EXPECT_DEATH(registry->transitionToCold(1), "Can only transition from HotElimination to ColdElimination");
    
    // Cannot mark ready for replacement from non-ColdElimination
    EXPECT_DEATH(registry->markReadyForReplacement(1), "Can only mark ColdElimination genomes as ready for replacement");
    
    // Cannot mark non-Active as Elite
    registry->markAsElite(2); // Active -> Elite (valid)
    EXPECT_DEATH(registry->markAsElite(2), "Can only mark Active genomes as Elite");
}

TEST_F(GlobalIndexRegistryTest, EliteProtection_StateTransitions) {
    // Test: Elite genomes have special transition rules
    
    // Elite can be marked for elimination
    registry->markAsElite(0);
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);
    registry->markForElimination(0);
    EXPECT_EQ(registry->getState(0), GenomeState::HotElimination);
    
    // Elite cannot be directly transitioned to other states
    registry->markAsElite(1);
    EXPECT_DEATH(registry->transitionToCold(1), "Can only transition from HotElimination to ColdElimination");
    EXPECT_DEATH(registry->markReadyForReplacement(1), "Can only mark ColdElimination genomes as ready for replacement");
}

// ============================================================================
// Index Lifecycle Management Tests
// ============================================================================

TEST_F(GlobalIndexRegistryTest, IndexReuse_LifecycleManagement) {
    // Test: getFreeIndex() properly reuses ReadyForReplacement indices
    
    // Initially no free indices
    EXPECT_EQ(registry->getFreeIndex(), INVALID_INDEX);
    
    // Create a ready-for-replacement index through full lifecycle
    registry->markForElimination(0);
    registry->transitionToCold(0);
    registry->markReadyForReplacement(0);
    EXPECT_EQ(registry->getState(0), GenomeState::ReadyForReplacement);
    
    // getFreeIndex should return this index and transition it to Active
    uint32_t freeIndex = registry->getFreeIndex();
    EXPECT_EQ(freeIndex, 0);
    EXPECT_EQ(registry->getState(0), GenomeState::Active);
    
    // No more free indices
    EXPECT_EQ(registry->getFreeIndex(), INVALID_INDEX);
}

TEST_F(GlobalIndexRegistryTest, IndexReuse_MultipleIndices) {
    // Test: Multiple ready-for-replacement indices are properly managed
    
    // Create multiple ready-for-replacement indices
    for (uint32_t i = 0; i < 3; ++i) {
        registry->markForElimination(i);
        registry->transitionToCold(i);
        registry->markReadyForReplacement(i);
    }
    
    // Get free indices in order
    std::vector<uint32_t> reusedIndices;
    uint32_t freeIndex;
    while ((freeIndex = registry->getFreeIndex()) != INVALID_INDEX) {
        reusedIndices.push_back(freeIndex);
        EXPECT_EQ(registry->getState(freeIndex), GenomeState::Active);
    }
    
    // Should have reused indices 0, 1, 2 in order
    EXPECT_EQ(reusedIndices.size(), 3u);
    EXPECT_EQ(reusedIndices[0], 0u);
    EXPECT_EQ(reusedIndices[1], 1u);
    EXPECT_EQ(reusedIndices[2], 2u);
}

TEST_F(GlobalIndexRegistryTest, EliteManagement_ProtectionAndClearing) {
    // Test: Elite status management and clearing
    
    // Mark multiple genomes as elite
    registry->markAsElite(0);
    registry->markAsElite(2);
    registry->markAsElite(4);
    
    // Verify elite status
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);
    EXPECT_EQ(registry->getState(1), GenomeState::Active);
    EXPECT_EQ(registry->getState(2), GenomeState::Elite);
    EXPECT_EQ(registry->getState(3), GenomeState::Active);
    EXPECT_EQ(registry->getState(4), GenomeState::Elite);
    
    // Clear all elite status
    registry->clearAllEliteStatus();
    
    // All elites should become active, others unchanged
    EXPECT_EQ(registry->getState(0), GenomeState::Active);
    EXPECT_EQ(registry->getState(1), GenomeState::Active);
    EXPECT_EQ(registry->getState(2), GenomeState::Active);
    EXPECT_EQ(registry->getState(3), GenomeState::Active);
    EXPECT_EQ(registry->getState(4), GenomeState::Active);
}

TEST_F(GlobalIndexRegistryTest, EliteManagement_PreservesOtherStates) {
    // Test: clearAllEliteStatus() doesn't affect non-elite states
    
    // Create mixed states
    registry->markAsElite(0);                    // Elite
    registry->markForElimination(1);             // HotElimination
    registry->markForElimination(2);
    registry->transitionToCold(2);               // ColdElimination
    registry->markForElimination(3);
    registry->transitionToCold(3);
    registry->markReadyForReplacement(3);        // ReadyForReplacement
    // Index 4 remains Active
    
    // Clear elite status
    registry->clearAllEliteStatus();
    
    // Only elite should change to active
    EXPECT_EQ(registry->getState(0), GenomeState::Active);           // Elite -> Active
    EXPECT_EQ(registry->getState(1), GenomeState::HotElimination);   // Unchanged
    EXPECT_EQ(registry->getState(2), GenomeState::ColdElimination);  // Unchanged
    EXPECT_EQ(registry->getState(3), GenomeState::ReadyForReplacement); // Unchanged
    EXPECT_EQ(registry->getState(4), GenomeState::Active);           // Unchanged
}

// ============================================================================
// Bounds Checking Tests
// ============================================================================

TEST_F(GlobalIndexRegistryTest, BoundsChecking_AssertionFailures) {
    // Test: All operations properly check index bounds
    
    uint32_t maxIndex = registry->getMaxIndex();
    EXPECT_EQ(maxIndex, 5u);
    
    // Out of bounds access should assert
    EXPECT_DEATH(registry->getState(maxIndex), "Global index out of range");
    EXPECT_DEATH(registry->markAsElite(maxIndex), "Global index out of range");
    EXPECT_DEATH(registry->markForElimination(maxIndex), "Global index out of range");
    EXPECT_DEATH(registry->transitionToCold(maxIndex), "Global index out of range");
    EXPECT_DEATH(registry->markReadyForReplacement(maxIndex), "Global index out of range");
}

TEST_F(GlobalIndexRegistryTest, BoundsChecking_ValidAccess) {
    // Test: All valid indices work correctly
    
    uint32_t maxIndex = registry->getMaxIndex();
    
    // All indices within bounds should work
    for (uint32_t i = 0; i < maxIndex; ++i) {
        EXPECT_EQ(registry->getState(i), GenomeState::Active);
        registry->markAsElite(i);
        EXPECT_EQ(registry->getState(i), GenomeState::Elite);
    }
}

// ============================================================================
// Registry Size and Bounds Tests
// ============================================================================

TEST_F(GlobalIndexRegistryTest, RegistrySize_InitialConfiguration) {
    // Test: Registry is properly initialized with correct size
    // Note: incrementMaxIndex() is protected and only accessible via PopulationContainer
    // This test focuses on the public interface and initial state
    
    uint32_t initialSize = registry->getMaxIndex();
    EXPECT_EQ(initialSize, 5u);
    
    // All initial indices should be Active and functional
    for (uint32_t i = 0; i < initialSize; ++i) {
        EXPECT_EQ(registry->getState(i), GenomeState::Active);
        registry->markAsElite(i);
        EXPECT_EQ(registry->getState(i), GenomeState::Elite);
    }
    
    // Clear all elites to restore Active state
    registry->clearAllEliteStatus();
    for (uint32_t i = 0; i < initialSize; ++i) {
        EXPECT_EQ(registry->getState(i), GenomeState::Active);
    }
}

// ============================================================================
// State Query Consistency Tests
// ============================================================================

TEST_F(GlobalIndexRegistryTest, StateQuery_Persistence) {
    // Test: getState() always returns accurate current state
    
    // Verify initial state
    EXPECT_EQ(registry->getState(0), GenomeState::Active);
    
    // State persists through multiple queries
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(registry->getState(0), GenomeState::Active);
    }
    
    // State changes are immediately reflected
    registry->markAsElite(0);
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);
    
    registry->markForElimination(0);
    EXPECT_EQ(registry->getState(0), GenomeState::HotElimination);
    
    registry->transitionToCold(0);
    EXPECT_EQ(registry->getState(0), GenomeState::ColdElimination);
    
    registry->markReadyForReplacement(0);
    EXPECT_EQ(registry->getState(0), GenomeState::ReadyForReplacement);
}

TEST_F(GlobalIndexRegistryTest, StateQuery_IndependentIndices) {
    // Test: State changes to one index don't affect others
    
    // Set different states for different indices
    registry->markAsElite(0);
    registry->markForElimination(1);
    registry->markForElimination(2);
    registry->transitionToCold(2);
    
    // Verify states are independent
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);
    EXPECT_EQ(registry->getState(1), GenomeState::HotElimination);
    EXPECT_EQ(registry->getState(2), GenomeState::ColdElimination);
    EXPECT_EQ(registry->getState(3), GenomeState::Active);
    EXPECT_EQ(registry->getState(4), GenomeState::Active);
    
    // Further changes maintain independence
    registry->transitionToCold(1);
    EXPECT_EQ(registry->getState(0), GenomeState::Elite);        // Unchanged
    EXPECT_EQ(registry->getState(1), GenomeState::ColdElimination); // Changed
    EXPECT_EQ(registry->getState(2), GenomeState::ColdElimination); // Unchanged
    EXPECT_EQ(registry->getState(3), GenomeState::Active);       // Unchanged
    EXPECT_EQ(registry->getState(4), GenomeState::Active);       // Unchanged
}