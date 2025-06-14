#include <gtest/gtest.h>
#include <vector>

#include "tests/test_common.h"
#include "version3/population/ResolveParentIndices.hpp"
#include "version3/population/ReproductiveInstruction.hpp"

using namespace Population;

class ResolveParentIndicesTest : public ::testing::Test {
protected:
    void SetUp() override {
        // No special setup needed for this simple operator
    }
    
    // Helper method to create instructions
    ReproductiveInstruction createInstruction(OperationType type, std::vector<uint32_t> relativeIndices) {
        ReproductiveInstruction instruction;
        instruction.operationType = type;
        instruction.relativeParentIndices = relativeIndices;
        instruction.globalParentIndices.clear(); // Start with empty global indices
        return instruction;
    }
    
    // Helper method to verify resolution
    void verifyResolution(const ReproductiveInstruction& instruction, 
                         const std::vector<uint32_t>& expectedGlobal) {
        EXPECT_EQ(instruction.globalParentIndices.size(), expectedGlobal.size()) 
            << "Global indices size mismatch";
        
        for (size_t i = 0; i < expectedGlobal.size(); ++i) {
            EXPECT_EQ(instruction.globalParentIndices[i], expectedGlobal[i])
                << "Global index mismatch at position " << i;
        }
    }
};

// =============================================================================
// CORE FUNCTIONALITY TESTS
// =============================================================================

TEST_F(ResolveParentIndicesTest, TestSingleInstructionSingleParent) {
    // Setup: 1 PRESERVE instruction with relativeParentIndices = [2]
    // GlobalIndices: [10, 20, 30, 40, 50]
    // Expected: globalParentIndices = [30]
    
    SpeciesInstructionSet instructions;
    instructions.push_back(createInstruction(OperationType::PRESERVE, {2}));
    
    std::vector<size_t> globalIndices = {10, 20, 30, 40, 50};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {30});
}

TEST_F(ResolveParentIndicesTest, TestSingleInstructionTwoParents) {
    // Setup: 1 CROSSOVER instruction with relativeParentIndices = [0, 3]
    // GlobalIndices: [10, 20, 30, 40, 50]
    // Expected: globalParentIndices = [10, 40]
    
    SpeciesInstructionSet instructions;
    instructions.push_back(createInstruction(OperationType::CROSSOVER, {0, 3}));
    
    std::vector<size_t> globalIndices = {10, 20, 30, 40, 50};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {10, 40});
}

TEST_F(ResolveParentIndicesTest, TestMultipleInstructionsIndependentProcessing) {
    // Setup: 3 instructions with different relativeParentIndices
    // Instruction 1: [1] → should resolve to [globalIndices[1]]
    // Instruction 2: [0, 2] → should resolve to [globalIndices[0], globalIndices[2]]
    // Instruction 3: [4] → should resolve to [globalIndices[4]]
    
    SpeciesInstructionSet instructions;
    instructions.push_back(createInstruction(OperationType::PRESERVE, {1}));
    instructions.push_back(createInstruction(OperationType::CROSSOVER, {0, 2}));
    instructions.push_back(createInstruction(OperationType::MUTATE_UNPROTECTED, {4}));
    
    std::vector<size_t> globalIndices = {100, 200, 300, 400, 500};
    
    resolveParentIndices(instructions, globalIndices);
    
    // Verify each instruction processed correctly and independently
    verifyResolution(instructions[0], {200});        // [1] → [200]
    verifyResolution(instructions[1], {100, 300});   // [0, 2] → [100, 300]
    verifyResolution(instructions[2], {500});        // [4] → [500]
}

// =============================================================================
// INSTRUCTION TYPE COVERAGE TESTS
// =============================================================================

TEST_F(ResolveParentIndicesTest, TestPreserveOperation) {
    // Setup: PRESERVE instruction with relativeParentIndices = [1]
    // Verify: Single parent index correctly resolved
    
    SpeciesInstructionSet instructions;
    instructions.push_back(ReproductiveInstruction::preserve(1));
    
    std::vector<size_t> globalIndices = {10, 20, 30};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {20});
}

TEST_F(ResolveParentIndicesTest, TestCrossoverOperation) {
    // Setup: CROSSOVER instruction with relativeParentIndices = [2, 4]
    // Verify: Both parent indices correctly resolved
    
    SpeciesInstructionSet instructions;
    instructions.push_back(ReproductiveInstruction::crossover(2, 4));
    
    std::vector<size_t> globalIndices = {10, 20, 30, 40, 50, 60};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {30, 50});
}

TEST_F(ResolveParentIndicesTest, TestMutateOperation) {
    // Setup: MUTATE_UNPROTECTED instruction with relativeParentIndices = [0]
    // Verify: Single parent index correctly resolved
    
    SpeciesInstructionSet instructions;
    instructions.push_back(ReproductiveInstruction::mutateUnprotected(0));
    
    std::vector<size_t> globalIndices = {15, 25, 35};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {15});
}

// =============================================================================
// EDGE CASE TESTS
// =============================================================================

TEST_F(ResolveParentIndicesTest, TestBoundaryIndices) {
    // Setup: Instructions with relative indices [0] and [size-1]
    // GlobalIndices: [100, 200, 300]
    // Test relative index 0 → global index 100
    // Test relative index 2 → global index 300
    
    SpeciesInstructionSet instructions;
    instructions.push_back(createInstruction(OperationType::PRESERVE, {0}));
    instructions.push_back(createInstruction(OperationType::PRESERVE, {2}));
    
    std::vector<size_t> globalIndices = {100, 200, 300};
    
    resolveParentIndices(instructions, globalIndices);
    
    verifyResolution(instructions[0], {100}); // First index
    verifyResolution(instructions[1], {300}); // Last index
}

TEST_F(ResolveParentIndicesTest, TestEmptyInstructionSet) {
    // Setup: Empty SpeciesInstructionSet
    // GlobalIndices: [10, 20, 30]
    // Expected: Function completes without error, no modifications
    
    SpeciesInstructionSet instructions; // Empty
    std::vector<size_t> globalIndices = {10, 20, 30};
    
    // Should not crash or throw
    EXPECT_NO_THROW(resolveParentIndices(instructions, globalIndices));
    
    // Instructions should remain empty
    EXPECT_TRUE(instructions.empty());
}

TEST_F(ResolveParentIndicesTest, TestClearsPreviousGlobalIndices) {
    // Setup: Instruction with pre-existing globalParentIndices = [999, 888]
    // Setup: relativeParentIndices = [1], globalIndices = [10, 20, 30]
    // Expected: globalParentIndices = [20] (old values cleared)
    
    SpeciesInstructionSet instructions;
    auto instruction = createInstruction(OperationType::PRESERVE, {1});
    
    // Pre-populate with old values
    instruction.globalParentIndices = {999, 888};
    instructions.push_back(instruction);
    
    std::vector<size_t> globalIndices = {10, 20, 30};
    
    resolveParentIndices(instructions, globalIndices);
    
    // Should have cleared old values and set new one
    verifyResolution(instructions[0], {20});
}

// =============================================================================
// DEBUG ASSERTION TESTS (Debug Build Only)
// =============================================================================

#ifdef DEBUG
TEST_F(ResolveParentIndicesTest, TestOutOfBoundsAssertion) {
    // Setup: Instruction with relativeParentIndices = [5]
    // GlobalIndices: [10, 20, 30] (size = 3, so index 5 is invalid)
    // Expected: Debug assertion should fire
    
    SpeciesInstructionSet instructions;
    instructions.push_back(createInstruction(OperationType::PRESERVE, {5})); // Out of bounds
    
    std::vector<size_t> globalIndices = {10, 20, 30}; // Size 3, index 5 invalid
    
    // Should trigger debug assertion
    EXPECT_DEATH(resolveParentIndices(instructions, globalIndices), "Relative index out of bounds");
}
#endif

// =============================================================================
// INTEGRATION-STYLE TESTS
// =============================================================================

TEST_F(ResolveParentIndicesTest, TestRealisticSmallSpecies) {
    // Setup: 3 instructions representing small species (2-3 genomes)
    // GlobalIndices: [15, 25] (species of size 2)
    // Mix of PRESERVE, CROSSOVER, MUTATE with various relative indices
    
    SpeciesInstructionSet instructions;
    instructions.push_back(ReproductiveInstruction::preserve(0));          // Preserve best
    instructions.push_back(ReproductiveInstruction::crossover(0, 1));      // Cross best with second
    instructions.push_back(ReproductiveInstruction::mutateUnprotected(1)); // Mutate second
    
    std::vector<size_t> globalIndices = {15, 25}; // Small species
    
    resolveParentIndices(instructions, globalIndices);
    
    // Verify all resolve correctly
    verifyResolution(instructions[0], {15});     // preserve(0) → [15]
    verifyResolution(instructions[1], {15, 25}); // crossover(0,1) → [15, 25]
    verifyResolution(instructions[2], {25});     // mutate(1) → [25]
}

TEST_F(ResolveParentIndicesTest, TestRealisticLargerSpecies) {
    // Setup: Multiple instructions for larger species (8-10 genomes)
    // GlobalIndices: [100, 110, 120, 130, 140, 150, 160, 170]
    // Various instructions with indices throughout the range
    
    SpeciesInstructionSet instructions;
    instructions.push_back(ReproductiveInstruction::preserve(0));          // Best
    instructions.push_back(ReproductiveInstruction::preserve(1));          // Second best
    instructions.push_back(ReproductiveInstruction::crossover(0, 2));      // Cross best with third
    instructions.push_back(ReproductiveInstruction::crossover(1, 3));      // Cross second with fourth
    instructions.push_back(ReproductiveInstruction::mutateUnprotected(4)); // Mutate middle
    instructions.push_back(ReproductiveInstruction::mutateProtected(7));   // Mutate worst
    
    std::vector<size_t> globalIndices = {100, 110, 120, 130, 140, 150, 160, 170};
    
    resolveParentIndices(instructions, globalIndices);
    
    // Verify correct resolution and no index confusion
    verifyResolution(instructions[0], {100});     // preserve(0) → [100]
    verifyResolution(instructions[1], {110});     // preserve(1) → [110]
    verifyResolution(instructions[2], {100, 120}); // crossover(0,2) → [100, 120]
    verifyResolution(instructions[3], {110, 130}); // crossover(1,3) → [110, 130]
    verifyResolution(instructions[4], {140});     // mutate(4) → [140]
    verifyResolution(instructions[5], {170});     // mutate(7) → [170]
}

// =============================================================================
// PERFORMANCE/EFFICIENCY TESTS
// =============================================================================

TEST_F(ResolveParentIndicesTest, TestNoObjectCreation) {
    // Setup: Instructions and verify only in-place modification occurs
    // Check that instruction objects themselves aren't copied/moved
    // Verify globalParentIndices vectors are reused, not recreated
    
    SpeciesInstructionSet instructions;
    auto originalInstruction = createInstruction(OperationType::CROSSOVER, {0, 2});
    instructions.push_back(originalInstruction);
    
    // Get pointer to the instruction for identity verification
    auto* instructionPtr = &instructions[0];
    auto* globalIndicesPtr = &instructions[0].globalParentIndices;
    
    std::vector<size_t> globalIndices = {10, 20, 30, 40, 50};
    
    resolveParentIndices(instructions, globalIndices);
    
    // Verify instruction object identity unchanged (no copying)
    EXPECT_EQ(instructionPtr, &instructions[0]) << "Instruction object should not be copied";
    
    // Verify globalParentIndices vector identity unchanged (in-place modification)
    EXPECT_EQ(globalIndicesPtr, &instructions[0].globalParentIndices) 
        << "globalParentIndices vector should be modified in-place";
    
    // Verify correct resolution occurred
    verifyResolution(instructions[0], {10, 30});
}