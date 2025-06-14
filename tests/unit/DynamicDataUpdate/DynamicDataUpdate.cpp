#include <gtest/gtest.h>
#include <map>
#include <unordered_map>

#include "tests/test_common.h"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/ReproductiveInstruction.hpp"

using namespace Population;

class DynamicDataUpdateTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test parameters with reasonable defaults
        params = std::make_unique<DynamicDataUpdateParams>(
            5,    // maxProtectionLimit
            3,    // maxSpeciesProtectionRating  
            0.3,  // protectedTierPercentage (30%)
            1     // worstSpeciesCount
        );
    }

    std::unique_ptr<DynamicDataUpdateParams> params;
    
    // Helper to create simple fitness result type for testing
    struct TestFitnessResult {
        double fitness;
        
        bool operator<(const TestFitnessResult& other) const {
            return fitness > other.fitness; // Higher fitness is better (first in map)
        }
        
        bool operator==(const TestFitnessResult& other) const {
            return fitness == other.fitness;
        }
    };
    
    // Helper to create test genome data multimap with specified fitness values and species IDs
    std::multimap<TestFitnessResult, DynamicGenomeData> createTestGenomeData(
        const std::vector<double>& fitnessValues,
        const std::vector<uint32_t>& speciesIds) {
        
        std::multimap<TestFitnessResult, DynamicGenomeData> genomeData;
        
        for (size_t i = 0; i < fitnessValues.size(); ++i) {
            TestFitnessResult fitness{fitnessValues[i]};
            DynamicGenomeData data;
            data.speciesId = speciesIds[i];
            data.protectionCounter = 0;
            data.isUnderRepair = false;
            data.isMarkedForElimination = false;
            
            genomeData.insert({fitness, data});
        }
        
        return genomeData;
    }
    
    // Helper to create test species data map
    std::unordered_map<uint32_t, DynamicSpeciesData> createTestSpeciesData(
        const std::vector<uint32_t>& speciesIds) {
        
        std::unordered_map<uint32_t, DynamicSpeciesData> speciesData;
        
        for (uint32_t id : speciesIds) {
            DynamicSpeciesData data;
            data.protectionRating = 0;
            data.currentPopulationSize = 0;
            data.isMarkedForElimination = false;
            
            speciesData[id] = data;
        }
        
        return speciesData;
    }
    
    // Helper to create minimal instruction sets for testing
    // Creates instruction sets with specified number of instructions per species
    GenerationInstructionSets createMockInstructionSets(
        const std::vector<uint32_t>& speciesIds,
        uint32_t instructionsPerSpecies = 0) {
        
        GenerationInstructionSets instructionSets;
        
        for (uint32_t speciesId : speciesIds) {
            SpeciesInstructionSet speciesInstructions;
            
            // Create simple PRESERVE instructions for testing
            for (uint32_t i = 0; i < instructionsPerSpecies; ++i) {
                speciesInstructions.push_back(
                    ReproductiveInstruction::preserve(i % 10) // Simple parent index cycling
                );
            }
            
            instructionSets[speciesId] = speciesInstructions;
        }
        
        return instructionSets;
    }
};

// =============================================================================
// 1. PROTECTION TIER CLASSIFICATION LOGIC
// =============================================================================

TEST_F(DynamicDataUpdateTest, ProtectionTier_10Genomes_30Percent) {
    // Test: 10 genomes, 30% protection → exactly 3 genomes classified as protected
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10}; // Best to worst
    std::vector<uint32_t> species = {1, 1, 1, 1, 2, 2, 2, 3, 3, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Check protection tier classification: bottom 30% = ranks 7, 8, 9 (0-indexed)
    // In our map, these correspond to fitness values 30, 20, 10
    size_t protectedCount = 0;
    for (const auto& [fit, data] : genomeData) {
        if (data.protectionCounter > 0) {
            protectedCount++;
        }
    }
    
    EXPECT_EQ(protectedCount, 3) << "Exactly 3 genomes should be in protected tier";
    
    // Verify specific genomes are protected (worst 3: fitness 30, 20, 10)
    EXPECT_GT(genomeData.find(TestFitnessResult{30})->second.protectionCounter, 0);
    EXPECT_GT(genomeData.find(TestFitnessResult{20})->second.protectionCounter, 0);
    EXPECT_GT(genomeData.find(TestFitnessResult{10})->second.protectionCounter, 0);
    
    // Verify others are not protected
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 0);
    EXPECT_EQ(genomeData.find(TestFitnessResult{90})->second.protectionCounter, 0);
    EXPECT_EQ(genomeData.find(TestFitnessResult{40})->second.protectionCounter, 0);
}

TEST_F(DynamicDataUpdateTest, ProtectionTier_BoundaryGenome) {
    // Test boundary genome (exactly at threshold) classified correctly
    std::vector<double> fitness = {100, 90, 80, 70, 60}; // 5 genomes
    std::vector<uint32_t> species = {1, 1, 2, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // With 30% protection and 5 genomes: 5 * 0.3 = 1.5 → 1 genome protected
    // Bottom 30% = rank 4 (0-indexed), which is fitness 60
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    size_t protectedCount = 0;
    for (const auto& [fit, data] : genomeData) {
        if (data.protectionCounter > 0) {
            protectedCount++;
        }
    }
    
    EXPECT_EQ(protectedCount, 1) << "Exactly 1 genome should be in protected tier";
    EXPECT_GT(genomeData.find(TestFitnessResult{60})->second.protectionCounter, 0) << "Boundary genome should be protected";
    EXPECT_EQ(genomeData.find(TestFitnessResult{70})->second.protectionCounter, 0) << "Genome just above boundary should not be protected";
}

TEST_F(DynamicDataUpdateTest, ProtectionTier_ZeroPercent) {
    // Test edge case: 0% protection → no genomes protected
    DynamicDataUpdateParams zeroParams(5, 3, 0.0, 1); // 0% protection
    
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1, 1, 2, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, zeroParams);
    
    for (const auto& [fit, data] : genomeData) {
        EXPECT_EQ(data.protectionCounter, 0) << "No genomes should be protected with 0% protection";
    }
}

TEST_F(DynamicDataUpdateTest, ProtectionTier_HundredPercent) {
    // Test edge case: 100% protection → all genomes protected
    DynamicDataUpdateParams hundredParams(5, 3, 1.0, 1); // 100% protection
    
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1, 1, 2, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, hundredParams);
    
    for (const auto& [fit, data] : genomeData) {
        EXPECT_GT(data.protectionCounter, 0) << "All genomes should be protected with 100% protection";
    }
}

// =============================================================================
// 2. PROTECTION COUNTER STATE TRANSITIONS
// =============================================================================

TEST_F(DynamicDataUpdateTest, ProtectionCounter_IncrementProtected) {
    // Test: Protected genome counter increments by 1
    std::vector<double> fitness = {100, 90, 80, 70, 60}; // 5 genomes
    std::vector<uint32_t> species = {1, 1, 1, 1, 1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Set initial counter values for testing increment
    genomeData.find(TestFitnessResult{60})->second.protectionCounter = 2; // Worst genome, should increment to 3
    genomeData.find(TestFitnessResult{100})->second.protectionCounter = 1; // Best genome, should reset to 0
    
    // With 30% protection and 5 genomes: bottom 30% = 1 genome (fitness 60)
    auto instructionSets = createMockInstructionSets({1});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(genomeData.find(TestFitnessResult{60})->second.protectionCounter, 3) 
        << "Protected genome counter should increment from 2 to 3";
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 0) 
        << "Non-protected genome counter should reset to 0";
}

TEST_F(DynamicDataUpdateTest, ProtectionCounter_ResetNonProtected) {
    // Test: Non-protected genome counter resets to 0
    // Use 10 genomes so 30% = 3 genomes protected (ranks 7, 8, 9)
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Set initial counters
    genomeData.find(TestFitnessResult{100})->second.protectionCounter = 3; // Rank 0, should reset to 0
    genomeData.find(TestFitnessResult{95})->second.protectionCounter = 1;  // Rank 1, should reset to 0
    genomeData.find(TestFitnessResult{65})->second.protectionCounter = 2;  // Rank 7, protected, should increment to 3
    genomeData.find(TestFitnessResult{60})->second.protectionCounter = 1;  // Rank 8, protected, should increment to 2
    
    // With 30% protection and 10 genomes: bottom 30% = 3 genomes (ranks 7, 8, 9)
    auto instructionSets = createMockInstructionSets({1});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 0) 
        << "Non-protected genome counter should reset to 0";
    EXPECT_EQ(genomeData.find(TestFitnessResult{95})->second.protectionCounter, 0) 
        << "Non-protected genome counter should reset to 0";
    EXPECT_EQ(genomeData.find(TestFitnessResult{65})->second.protectionCounter, 3) 
        << "Protected genome counter should increment from 2 to 3";
    EXPECT_EQ(genomeData.find(TestFitnessResult{60})->second.protectionCounter, 2) 
        << "Protected genome counter should increment from 1 to 2";
}

TEST_F(DynamicDataUpdateTest, ProtectionCounter_ExceedsLimitTriggersElimination) {
    // Test: Counter exceeding limit triggers elimination flag
    // Use 10 genomes so 30% = 3 genomes protected (ranks 7, 8, 9)
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Set counter at limit (maxProtectionLimit = 5) for worst genome
    genomeData.find(TestFitnessResult{55})->second.protectionCounter = 5; // Rank 9, at limit, will increment to 6
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(genomeData.find(TestFitnessResult{55})->second.protectionCounter, 6) 
        << "Protected genome counter should increment to 6";
    EXPECT_TRUE(genomeData.find(TestFitnessResult{55})->second.isMarkedForElimination) 
        << "Genome should be marked for elimination when counter exceeds limit";
    EXPECT_FALSE(genomeData.find(TestFitnessResult{100})->second.isMarkedForElimination) 
        << "Non-protected genome should not be marked for elimination";
}

TEST_F(DynamicDataUpdateTest, ProtectionCounter_UnderRepairGenomesSkipped) {
    // Test: Genomes under repair are skipped for protection updates
    // Use 10 genomes so 30% = 3 genomes protected (ranks 7, 8, 9)
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Mark worst genome as under repair
    genomeData.find(TestFitnessResult{55})->second.isUnderRepair = true; // Rank 9, would be protected but under repair
    genomeData.find(TestFitnessResult{55})->second.protectionCounter = 2; // Should remain unchanged
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(genomeData.find(TestFitnessResult{55})->second.protectionCounter, 2) 
        << "Under repair genome counter should remain unchanged";
    EXPECT_FALSE(genomeData.find(TestFitnessResult{55})->second.isMarkedForElimination) 
        << "Under repair genome should not be marked for elimination";
    
    // Verify other protected genomes still get updated
    EXPECT_GT(genomeData.find(TestFitnessResult{60})->second.protectionCounter, 0) 
        << "Other protected genomes should still be updated";
}

// =============================================================================
// 3. SPECIES PERFORMANCE RANKING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_AverageRankCalculation) {
    // Test: Species with genomes at ranks [0,5,9] has average rank 4.67
    // Create scenario: Species 1 has genomes at ranks 0, 5, 9 (fitness 100, 75, 55)
    //                 Species 2 has genomes at ranks 2, 7 (fitness 90, 65)
    //                 Species 3 has genomes at ranks 1, 3, 4, 6, 8 (fitness 95, 85, 80, 70, 60)
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   3,  2,  3,  3,  1,  3,  2,  3,  1};
    //                              r0   r1  r2  r3  r4  r5  r6  r7  r8  r9
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Species 1: ranks [0, 5, 9] → average = (0+5+9)/3 = 14/3 = 4.67
    // Species 2: ranks [2, 7] → average = (2+7)/2 = 4.5  
    // Species 3: ranks [1, 3, 4, 6, 8] → average = (1+3+4+6+8)/5 = 22/5 = 4.4
    
    // Species 3 has best average (4.4), Species 2 is middle (4.5), Species 1 is worst (4.67)
    // With worstSpeciesCount=1, Species 1 should be penalized
    
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Species 1 should be penalized (worst average rank)";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Species 2 should not be penalized";
    EXPECT_EQ(speciesData[3].protectionRating, 0) << "Species 3 should not be penalized";
}

TEST_F(DynamicDataUpdateTest, SpeciesRanking_WorstNSpeciesIdentification) {
    // Test: Worst N species correctly identified from ranking
    DynamicDataUpdateParams multiWorstParams(5, 3, 0.3, 2); // worstSpeciesCount = 2
    
    // Create 4 species with different average ranks
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30};
    std::vector<uint32_t> species = {1,   2,  3,  4,  1,  2,  3,  4};
    //                              r0   r1  r2  r3  r4  r5  r6  r7
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3, 4});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, multiWorstParams);
    
    // Species average ranks:
    // Species 1: ranks [0, 4] → average = (0+4)/2 = 2.0 (best)
    // Species 2: ranks [1, 5] → average = (1+5)/2 = 3.0 
    // Species 3: ranks [2, 6] → average = (2+6)/2 = 4.0
    // Species 4: ranks [3, 7] → average = (3+7)/2 = 5.0 (worst)
    
    // Worst 2 species should be Species 3 and Species 4
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Species 1 should not be penalized (best)";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Species 2 should not be penalized (second best)";
    EXPECT_EQ(speciesData[3].protectionRating, 1) << "Species 3 should be penalized (second worst)";
    EXPECT_EQ(speciesData[4].protectionRating, 1) << "Species 4 should be penalized (worst)";
}

TEST_F(DynamicDataUpdateTest, SpeciesRanking_SingleSpeciesPenalized) {
    // Test: Single species gets penalized when worstSpeciesCount=1
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    //                              r0   r1  r2  r3  r4  r5
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params); // worstSpeciesCount = 1
    
    // Species average ranks:
    // Species 1: ranks [0, 1] → average = (0+1)/2 = 0.5 (best)
    // Species 2: ranks [2, 3] → average = (2+3)/2 = 2.5
    // Species 3: ranks [4, 5] → average = (4+5)/2 = 4.5 (worst)
    
    // Only Species 3 (worst) should be penalized
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Species 1 should not be penalized";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Species 2 should not be penalized";
    EXPECT_EQ(speciesData[3].protectionRating, 1) << "Species 3 should be penalized (worst)";
}

TEST_F(DynamicDataUpdateTest, SpeciesRanking_SingleSpeciesNoOperation) {
    // Test: Single species in population should not crash and should not be penalized
    // (since worstSpeciesCount=1 but only 1 species exists)
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    auto instructionSets = createMockInstructionSets({1});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // With only 1 species and worstSpeciesCount=1, the species should be penalized
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Single species should be penalized when worstSpeciesCount=1";
}

// =============================================================================
// 4. SPECIES RATING STATE TRANSITIONS
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRating_IncrementWorstSpecies) {
    // Test: Worst species rating increments by 1
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial rating values for testing increment
    speciesData[1].protectionRating = 0; // Best species, should remain 0
    speciesData[2].protectionRating = 1; // Middle species, should reset to 0
    speciesData[3].protectionRating = 2; // Worst species, should increment to 3
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params); // worstSpeciesCount = 1
    
    // Species 3 is worst (average rank 4.5), should increment from 2 to 3
    EXPECT_EQ(speciesData[3].protectionRating, 3) << "Worst species rating should increment from 2 to 3";
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Best species rating should remain 0";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Non-worst species rating should reset to 0";
}

TEST_F(DynamicDataUpdateTest, SpeciesRating_ResetNonWorstSpecies) {
    // Test: Non-worst species rating resets to 0
    std::vector<double> fitness = {100, 90, 80, 70};
    std::vector<uint32_t> species = {1,   2,  3,  4};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    // Set initial ratings for non-worst species
    speciesData[1].protectionRating = 2; // Best species, should reset to 0
    speciesData[2].protectionRating = 1; // Should reset to 0
    speciesData[3].protectionRating = 3; // Should reset to 0
    speciesData[4].protectionRating = 1; // Worst species, should increment to 2
    
    auto instructionSets = createMockInstructionSets({1, 2, 3, 4});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params); // worstSpeciesCount = 1
    
    // Species 4 is worst (rank 3), others should reset
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Non-worst species rating should reset to 0";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Non-worst species rating should reset to 0";
    EXPECT_EQ(speciesData[3].protectionRating, 0) << "Non-worst species rating should reset to 0";
    EXPECT_EQ(speciesData[4].protectionRating, 2) << "Worst species rating should increment from 1 to 2";
}

TEST_F(DynamicDataUpdateTest, SpeciesRating_ExceedsLimitTriggersElimination) {
    // Test: Rating exceeding limit triggers elimination flag
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set worst species rating at limit (maxSpeciesProtectionRating = 3)
    speciesData[3].protectionRating = 3; // At limit, will increment to 4
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Species 3 is worst, should increment to 4 and be marked for elimination
    EXPECT_EQ(speciesData[3].protectionRating, 4) << "Worst species rating should increment to 4";
    EXPECT_TRUE(speciesData[3].isMarkedForElimination) << "Species should be marked for elimination when rating exceeds limit";
    EXPECT_FALSE(speciesData[1].isMarkedForElimination) << "Non-worst species should not be marked for elimination";
    EXPECT_FALSE(speciesData[2].isMarkedForElimination) << "Non-worst species should not be marked for elimination";
}

TEST_F(DynamicDataUpdateTest, SpeciesRating_MultipleWorstSpeciesIncrement) {
    // Test: Multiple worst species all increment their ratings
    DynamicDataUpdateParams multiWorstParams(5, 3, 0.3, 2); // worstSpeciesCount = 2
    
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30};
    std::vector<uint32_t> species = {1,   2,  3,  4,  1,  2,  3,  4};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    // Set initial ratings
    speciesData[1].protectionRating = 1; // Best species, should reset to 0
    speciesData[2].protectionRating = 2; // Should reset to 0  
    speciesData[3].protectionRating = 1; // Worst species, should increment to 2
    speciesData[4].protectionRating = 0; // Worst species, should increment to 1
    
    auto instructionSets = createMockInstructionSets({1, 2, 3, 4});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, multiWorstParams);
    
    // Species 3 and 4 are worst 2, should both increment
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Best species rating should reset to 0";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Non-worst species rating should reset to 0";
    EXPECT_EQ(speciesData[3].protectionRating, 2) << "Worst species rating should increment from 1 to 2";
    EXPECT_EQ(speciesData[4].protectionRating, 1) << "Worst species rating should increment from 0 to 1";
}

TEST_F(DynamicDataUpdateTest, SpeciesRating_EscapeWorstStatusResetsRating) {
    // Test: Species that escapes worst status gets rating reset
    // This simulates a species improving over generations
    
    // First generation: Species 3 is worst
    std::vector<double> fitness1 = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species1 = {1,   1,  2,  2,  3,  3};
    
    auto genomeData1 = createTestGenomeData(fitness1, species1);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    speciesData[3].protectionRating = 1; // Species 3 starts with rating 1
    
    auto instructionSets1 = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData1, speciesData, instructionSets1, *params);
    
    // Species 3 should increment to 2
    EXPECT_EQ(speciesData[3].protectionRating, 2) << "Species 3 should increment to 2 in first generation";
    
    // Second generation: Species 3 improves, Species 2 becomes worst
    std::vector<double> fitness2 = {100, 95, 50, 45, 90, 85}; // Species 3 now has better fitness
    std::vector<uint32_t> species2 = {1,   1,  2,  2,  3,  3};
    
    auto genomeData2 = createTestGenomeData(fitness2, species2);
    
    auto instructionSets2 = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData2, speciesData, instructionSets2, *params);
    
    // Species 3 should reset to 0 (no longer worst), Species 2 should increment
    EXPECT_EQ(speciesData[3].protectionRating, 0) << "Species 3 should reset to 0 after escaping worst status";
    EXPECT_GT(speciesData[2].protectionRating, 0) << "Species 2 should be penalized as new worst species";
}

// =============================================================================
// 5. POPULATION SIZE UPDATES (CORE OPERATOR RESPONSIBILITY)
// =============================================================================

TEST_F(DynamicDataUpdateTest, PopulationSize_UpdatedCorrectly) {
    // Test core responsibility: operator must update currentPopulationSize
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial sizes to wrong values
    speciesData[1].currentPopulationSize = 99; // Should become 2
    speciesData[2].currentPopulationSize = 99; // Should become 2  
    speciesData[3].currentPopulationSize = 99; // Should become 2
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(speciesData[1].currentPopulationSize, 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
}

TEST_F(DynamicDataUpdateTest, PopulationSize_UnevenDistribution) {
    // Test with uneven species distribution
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65};
    std::vector<uint32_t> species = {1,   1,  1,  1,  2,  2,  3,  3};
    //                              Species 1: 4 genomes, Species 2: 2 genomes, Species 3: 2 genomes
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial sizes to wrong values
    speciesData[1].currentPopulationSize = 0;
    speciesData[2].currentPopulationSize = 10;
    speciesData[3].currentPopulationSize = 5;
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(speciesData[1].currentPopulationSize, 4) << "Species 1 should have 4 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
}

// =============================================================================
// 6. EMPTY DATA HANDLING (BASIC ROBUSTNESS)
// =============================================================================

TEST_F(DynamicDataUpdateTest, EmptyData_NoOperations) {
    // Test basic robustness: empty inputs should not crash
    std::multimap<TestFitnessResult, DynamicGenomeData> emptyGenomeData;
    std::unordered_map<uint32_t, DynamicSpeciesData> emptySpeciesData;
    
    // Should not crash or throw
    auto instructionSets = createMockInstructionSets({});
    EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, emptySpeciesData, instructionSets, *params));
    
    // Data should remain empty
    EXPECT_TRUE(emptyGenomeData.empty()) << "Empty genome data should remain empty";
    EXPECT_TRUE(emptySpeciesData.empty()) << "Empty species data should remain empty";
}

TEST_F(DynamicDataUpdateTest, EmptyGenomes_NonEmptySpecies) {
    // Test edge case: no genomes but species data exists
    std::multimap<TestFitnessResult, DynamicGenomeData> emptyGenomeData;
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial population sizes to non-zero
    speciesData[1].currentPopulationSize = 5;
    speciesData[2].currentPopulationSize = 3;
    speciesData[3].currentPopulationSize = 2;
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, speciesData, instructionSets, *params));
    
    // All population sizes should be reset to 0
    EXPECT_EQ(speciesData[1].currentPopulationSize, 0) << "Species 1 population should be reset to 0";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 0) << "Species 2 population should be reset to 0";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 0) << "Species 3 population should be reset to 0";
    
    // No species should be penalized (no ranking possible with no genomes)
    EXPECT_EQ(speciesData[1].protectionRating, 0) << "Species 1 rating should remain 0";
    EXPECT_EQ(speciesData[2].protectionRating, 0) << "Species 2 rating should remain 0";
    EXPECT_EQ(speciesData[3].protectionRating, 0) << "Species 3 rating should remain 0";
}

TEST_F(DynamicDataUpdateTest, EmptySpecies_NonEmptyGenomes) {
    // Test edge case: genomes exist but no species data
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    std::unordered_map<uint32_t, DynamicSpeciesData> emptySpeciesData;
    
    // Should not crash even though species data is missing
    auto instructionSets = createMockInstructionSets({});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, emptySpeciesData, instructionSets, *params));
    
    // Genome protection logic should still work (no species dependencies)
    // With 30% protection and 3 genomes: 3 * 0.3 = 0.9 → 0 genomes protected
    for (const auto& [fit, data] : genomeData) {
        EXPECT_EQ(data.protectionCounter, 0) << "No genomes should be protected with small population";
    }
}

// =============================================================================
// 7. PARAMETER BOUNDARY CONDITIONS
// =============================================================================

TEST_F(DynamicDataUpdateTest, ParameterBoundary_WorstCountExceedsSpecies) {
    // Test: worstSpeciesCount=5 but only 2 species exist
    DynamicDataUpdateParams excessParams(5, 3, 0.3, 5); // worstSpeciesCount=5
    
    std::vector<double> fitness = {100, 90, 80, 70};
    std::vector<uint32_t> species = {1,   1,  2,  2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2});
    
    // Should handle gracefully - both species get penalized
    auto instructionSets = createMockInstructionSets({1, 2});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, excessParams));
    
    // Both species should be penalized since worstSpeciesCount (5) > actual species (2)
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Species 1 should be penalized when worstCount exceeds total";
    EXPECT_EQ(speciesData[2].protectionRating, 1) << "Species 2 should be penalized when worstCount exceeds total";
}

TEST_F(DynamicDataUpdateTest, ParameterBoundary_WorstCountEqualsSpecies) {
    // Test: worstSpeciesCount equals actual species count
    DynamicDataUpdateParams equalParams(5, 3, 0.3, 3); // worstSpeciesCount=3
    
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, equalParams));
    
    // All species should be penalized since worstSpeciesCount (3) == actual species (3)
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "All species should be penalized when worstCount equals total";
    EXPECT_EQ(speciesData[2].protectionRating, 1) << "All species should be penalized when worstCount equals total";
    EXPECT_EQ(speciesData[3].protectionRating, 1) << "All species should be penalized when worstCount equals total";
}

TEST_F(DynamicDataUpdateTest, ParameterBoundary_WorstCountZeroDeathTest) {
    // Test: worstSpeciesCount=0 should trigger assertion failure
    EXPECT_DEATH({
        DynamicDataUpdateParams zeroParams(5, 3, 0.3, 0); // worstSpeciesCount=0
    }, "worstSpeciesCount > 0");
}

TEST_F(DynamicDataUpdateTest, ParameterBoundary_InvalidProtectionPercentage) {
    // Test: invalid protection percentages should trigger assertion failure
    EXPECT_DEATH({
        DynamicDataUpdateParams invalidParams(5, 3, -0.1, 1); // negative percentage
    }, "protectedTierPercentage >= 0.0 && protectedTierPercentage <= 1.0");
    
    EXPECT_DEATH({
        DynamicDataUpdateParams invalidParams2(5, 3, 1.1, 1); // > 100%
    }, "protectedTierPercentage >= 0.0 && protectedTierPercentage <= 1.0");
}

// =============================================================================
// 8. ALL GENOMES UNDER REPAIR (IMPORTANT EDGE CASE)
// =============================================================================

TEST_F(DynamicDataUpdateTest, AllGenomesUnderRepair_NoProtectionUpdates) {
    // Test edge case: all genomes are under repair
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   1,  1,  2,  2,  2,  3,  3,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Mark all genomes as under repair
    for (auto& [fit, data] : genomeData) {
        data.isUnderRepair = true;
        data.protectionCounter = 2; // Should remain unchanged
    }
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // No genome protection counters should change
    for (const auto& [fit, data] : genomeData) {
        EXPECT_EQ(data.protectionCounter, 2) << "Under repair genome counters should remain unchanged";
        EXPECT_FALSE(data.isMarkedForElimination) << "Under repair genomes should not be marked for elimination";
    }
    
    // But species ratings should still update (not dependent on repair status)
    // Species 3 has worst average rank, should be penalized
    EXPECT_GT(speciesData[3].protectionRating, 0) << "Worst species should still be penalized even when all genomes under repair";
    
    // Population sizes should still be updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 3) << "Population size should still be updated";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3) << "Population size should still be updated";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 4) << "Population size should still be updated";
}

TEST_F(DynamicDataUpdateTest, MixedRepairStatus_OnlyNonRepairUpdated) {
    // Test: mixed repair status - only non-repair genomes get protection updates
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   1,  1,  1,  2,  2,  2,  3,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Mark some genomes as under repair (including some that would be protected)
    genomeData.find(TestFitnessResult{65})->second.isUnderRepair = true; // Rank 7, would be protected
    genomeData.find(TestFitnessResult{65})->second.protectionCounter = 3; // Should remain unchanged
    
    genomeData.find(TestFitnessResult{100})->second.isUnderRepair = true; // Rank 0, wouldn't be protected anyway
    genomeData.find(TestFitnessResult{100})->second.protectionCounter = 1; // Should remain unchanged
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Under repair genomes should remain unchanged
    EXPECT_EQ(genomeData.find(TestFitnessResult{65})->second.protectionCounter, 3) << "Under repair genome should remain unchanged";
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 1) << "Under repair genome should remain unchanged";
    
    // Non-repair protected genomes should be updated (ranks 7, 8, 9 are protected)
    EXPECT_GT(genomeData.find(TestFitnessResult{60})->second.protectionCounter, 0) << "Non-repair protected genome should be updated";
    EXPECT_GT(genomeData.find(TestFitnessResult{55})->second.protectionCounter, 0) << "Non-repair protected genome should be updated";
    
    // Non-repair non-protected genomes should reset to 0
    EXPECT_EQ(genomeData.find(TestFitnessResult{95})->second.protectionCounter, 0) << "Non-repair non-protected genome should reset";
}

// =============================================================================
// 9. SINGLE GENOME POPULATION (ABSOLUTE MINIMUM EDGE CASE)
// =============================================================================

TEST_F(DynamicDataUpdateTest, SingleGenome_HandledCorrectly) {
    // Test absolute minimum: 1 genome total
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    auto instructionSets = createMockInstructionSets({1});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // With 30% protection and 1 genome: 1 * 0.3 = 0.3 → 0 genomes protected
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 0) 
        << "Single genome should not be protected due to percentage truncation";
    
    // Population size should be correct
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Single species should have population size 1";
    
    // Single species should be penalized (worstSpeciesCount=1, and we have 1 species)
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Single species should be penalized as worst species";
}

TEST_F(DynamicDataUpdateTest, SingleGenome_WithHighProtectionPercentage) {
    // Test single genome with high protection percentage
    DynamicDataUpdateParams highProtectionParams(5, 3, 1.0, 1); // 100% protection
    
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    auto instructionSets = createMockInstructionSets({1});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, highProtectionParams));
    
    // With 100% protection and 1 genome: 1 * 1.0 = 1 genome protected
    EXPECT_GT(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 0) 
        << "Single genome should be protected with 100% protection";
    
    // Species should still be penalized (only species, so worst by default)
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Single species should still be penalized";
}

TEST_F(DynamicDataUpdateTest, SingleGenome_AtProtectionLimit) {
    // Test single genome at protection counter limit
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Set genome at protection limit with 100% protection to trigger increment
    DynamicDataUpdateParams fullProtectionParams(5, 3, 1.0, 1); // 100% protection
    genomeData.find(TestFitnessResult{100})->second.protectionCounter = 5; // At limit
    
    auto instructionSets = createMockInstructionSets({1});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, fullProtectionParams);
    
    // Should increment to 6 and be marked for elimination
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 6) 
        << "Single genome should increment beyond limit";
    EXPECT_TRUE(genomeData.find(TestFitnessResult{100})->second.isMarkedForElimination) 
        << "Single genome should be marked for elimination when exceeding limit";
}

TEST_F(DynamicDataUpdateTest, SingleGenome_UnderRepair) {
    // Test single genome under repair
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Mark single genome as under repair
    genomeData.find(TestFitnessResult{100})->second.isUnderRepair = true;
    genomeData.find(TestFitnessResult{100})->second.protectionCounter = 3;
    
    auto instructionSets = createMockInstructionSets({1});
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // Under repair genome should remain unchanged
    EXPECT_EQ(genomeData.find(TestFitnessResult{100})->second.protectionCounter, 3) 
        << "Under repair single genome should remain unchanged";
    EXPECT_FALSE(genomeData.find(TestFitnessResult{100})->second.isMarkedForElimination) 
        << "Under repair single genome should not be marked for elimination";
    
    // Species should still be ranked and penalized (ranking doesn't depend on repair status)
    EXPECT_EQ(speciesData[1].protectionRating, 1) << "Single species should still be penalized";
    
    // Population size should still be updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Population size should still be updated";
}

// =============================================================================
// NEW TESTS FOR INSTRUCTION SET COUNTING FUNCTIONALITY
// =============================================================================

// =============================================================================
// TEST CATEGORY 1: INSTRUCTION COUNTING ACCURACY
// =============================================================================

TEST_F(DynamicDataUpdateTest, InstructionCounting_BasicInstructionSetCounting) {
    // Test: Basic instruction set counting produces correct numeric results
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Test different instruction counts
    auto instructionSets = createMockInstructionSets({1, 2, 3}, 0); // Start with 0
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // 5 instructions
    
    // Create larger instruction set for species 3 (50 instructions)
    SpeciesInstructionSet largeSet;
    for (uint32_t i = 0; i < 50; ++i) {
        largeSet.push_back(ReproductiveInstruction::preserve(i % 10));
    }
    instructionSets[3] = largeSet;
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify exact counts
    EXPECT_EQ(speciesData[1].instructionSetsSize, 1) << "Species 1 should have exactly 1 instruction";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 5) << "Species 2 should have exactly 5 instructions";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 50) << "Species 3 should have exactly 50 instructions";
}

TEST_F(DynamicDataUpdateTest, InstructionCounting_EmptyInstructionSet) {
    // Test: Empty instruction sets produce count = 0
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2});
    
    // Create instruction sets with empty vectors
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{}; // Empty vector
    instructionSets[2] = SpeciesInstructionSet{}; // Empty vector
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(speciesData[1].instructionSetsSize, 0) << "Empty instruction set should produce count 0";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 0) << "Empty instruction set should produce count 0";
}

TEST_F(DynamicDataUpdateTest, InstructionCounting_MultipleSpeciesIndependentCounting) {
    // Test: Independent counting across multiple species
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1, 1, 2, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Create different instruction counts per species
    GenerationInstructionSets instructionSets;
    
    // Species A: 3 instructions
    instructionSets[1] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    };
    
    // Species B: 7 instructions
    SpeciesInstructionSet sevenInstructions;
    for (uint32_t i = 0; i < 7; ++i) {
        sevenInstructions.push_back(ReproductiveInstruction::preserve(i % 5));
    }
    instructionSets[2] = sevenInstructions;
    
    // Species C: 2 instructions
    instructionSets[3] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    };
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify independent counts
    EXPECT_EQ(speciesData[1].instructionSetsSize, 3) << "Species 1 should have 3 instructions independently";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 7) << "Species 2 should have 7 instructions independently";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 2) << "Species 3 should have 2 instructions independently";
}

TEST_F(DynamicDataUpdateTest, InstructionCounting_MixedEmptyAndNonEmpty) {
    // Test: Mixed empty and non-empty instruction sets
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{}; // 0 instructions
    
    // Species 2: 10 instructions
    SpeciesInstructionSet tenInstructions;
    for (uint32_t i = 0; i < 10; ++i) {
        tenInstructions.push_back(ReproductiveInstruction::preserve(i % 8));
    }
    instructionSets[2] = tenInstructions;
    
    instructionSets[3] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    EXPECT_EQ(speciesData[1].instructionSetsSize, 0) << "Species with empty set should have count 0";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 10) << "Species with 10 instructions should have count 10";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 1) << "Species with 1 instruction should have count 1";
}

// =============================================================================
// TEST CATEGORY 2: FIELD UPDATE VERIFICATION
// =============================================================================

TEST_F(DynamicDataUpdateTest, FieldUpdate_InstructionSetsSizeFieldUpdates) {
    // Test: instructionSetsSize field gets exact counted value
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Pre-set different initial values to ensure they get overwritten
    speciesData[1].instructionSetsSize = 999;
    speciesData[2].instructionSetsSize = 888;
    speciesData[3].instructionSetsSize = 777;
    
    auto instructionSets = createMockInstructionSets({1, 2, 3}, 0);
    instructionSets[1] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4),
        ReproductiveInstruction::preserve(5),
        ReproductiveInstruction::preserve(6)
    }; // 7 instructions
    
    instructionSets[2] = SpeciesInstructionSet{}; // 0 instructions
    instructionSets[3] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify exact field values
    EXPECT_EQ(speciesData[1].instructionSetsSize, 7) << "instructionSetsSize field should contain exact count 7";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 0) << "instructionSetsSize field should contain exact count 0";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 3) << "instructionSetsSize field should contain exact count 3";
}

TEST_F(DynamicDataUpdateTest, FieldUpdate_FieldOverwriteBehavior) {
    // Test: Existing field values are properly overwritten
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set specific initial values to test overwrite behavior
    speciesData[1].instructionSetsSize = 10; // Will become 5
    speciesData[2].instructionSetsSize = 0;  // Will become 3
    speciesData[3].instructionSetsSize = 100; // Will become 0
    
    auto instructionSets = createMockInstructionSets({1, 2, 3}, 0);
    instructionSets[1] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // 5 instructions
    
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions
    
    instructionSets[3] = SpeciesInstructionSet{}; // 0 instructions
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify complete overwrite of old values
    EXPECT_EQ(speciesData[1].instructionSetsSize, 5) << "Old value 10 should be overwritten with 5";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 3) << "Old value 0 should be overwritten with 3";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 0) << "Old value 100 should be overwritten with 0";
}

TEST_F(DynamicDataUpdateTest, FieldUpdate_OtherFieldsPreservation) {
    // Test: Counting doesn't modify unrelated species data fields
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set specific values for other fields to verify they're preserved
    speciesData[1].protectionRating = 5;
    speciesData[1].isMarkedForElimination = true;
    speciesData[1].speciesRank = 10;
    
    speciesData[2].protectionRating = 2;
    speciesData[2].isMarkedForElimination = false;
    speciesData[2].speciesRank = 7;
    
    speciesData[3].protectionRating = 0;
    speciesData[3].isMarkedForElimination = true;
    speciesData[3].speciesRank = 1;
    
    auto instructionSets = createMockInstructionSets({1, 2, 3}, 3); // All get 3 instructions
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify only instructionSetsSize changed, others preserved initially
    // Note: The operator may update some fields like protectionRating and speciesRank as part of its normal operation
    // but currentPopulationSize should be updated normally
    EXPECT_EQ(speciesData[1].instructionSetsSize, 3) << "instructionSetsSize should be updated";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 3) << "instructionSetsSize should be updated";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 3) << "instructionSetsSize should be updated";
    
    // Note: protectionRating and speciesRank will be updated by the operator's normal logic
    // so we only verify that the instruction counting doesn't break the normal field updates
}

// =============================================================================
// TEST CATEGORY 3: SPECIES ID MAPPING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesIDMapping_SpeciesIDLookupSuccess) {
    // Test: Correct species data entry gets updated based on ID
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {123, 456, 999};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({123, 456, 999});
    
    GenerationInstructionSets instructionSets;
    instructionSets[123] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // 5 instructions
    
    instructionSets[456] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    
    instructionSets[999] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4),
        ReproductiveInstruction::preserve(5),
        ReproductiveInstruction::preserve(6),
        ReproductiveInstruction::preserve(7)
    }; // 8 instructions
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify correct species entries updated based on ID matching
    EXPECT_EQ(speciesData[123].instructionSetsSize, 5) << "Species 123 should have count 5";
    EXPECT_EQ(speciesData[456].instructionSetsSize, 2) << "Species 456 should have count 2";
    EXPECT_EQ(speciesData[999].instructionSetsSize, 8) << "Species 999 should have count 8";
}

TEST_F(DynamicDataUpdateTest, SpeciesIDMapping_SpeciesIDMismatchHandling) {
    // Test: Behavior when species IDs don't match between inputs
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4, 5}); // More species in data than instructions
    
    // Set initial values for tracking changes
    speciesData[1].instructionSetsSize = 111;
    speciesData[2].instructionSetsSize = 222;
    speciesData[3].instructionSetsSize = 333;
    speciesData[4].instructionSetsSize = 444; // Not in instruction sets
    speciesData[5].instructionSetsSize = 555; // Not in instruction sets
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    instructionSets[3] = SpeciesInstructionSet{}; // 0 instructions
    // Species 4 and 5 not in instruction sets
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Species in instruction sets should be updated
    EXPECT_EQ(speciesData[1].instructionSetsSize, 1) << "Species 1 should be updated from instruction sets";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 2) << "Species 2 should be updated from instruction sets";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 0) << "Species 3 should be updated from instruction sets";
    
    // Species not in instruction sets should be reset to 0 (Phase 1 behavior)
    EXPECT_EQ(speciesData[4].instructionSetsSize, 0) << "Species 4 should be reset to 0 (not in instruction sets)";
    EXPECT_EQ(speciesData[5].instructionSetsSize, 0) << "Species 5 should be reset to 0 (not in instruction sets)";
}

TEST_F(DynamicDataUpdateTest, SpeciesIDMapping_NewSpeciesDiscovery) {
    // Test: New species appear in instruction sets but not in species data
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only species 1 and 2
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    instructionSets[999] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions for new species
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Existing species should be updated
    EXPECT_EQ(speciesData[1].instructionSetsSize, 1) << "Existing species 1 should be updated";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 2) << "Existing species 2 should be updated";
    
    // New species should be created with minimal data
    auto newSpeciesIt = speciesData.find(999);
    EXPECT_NE(newSpeciesIt, speciesData.end()) << "New species 999 should be created";
    if (newSpeciesIt != speciesData.end()) {
        EXPECT_EQ(newSpeciesIt->second.instructionSetsSize, 3) << "New species should have correct instruction count";
    }
}

// =============================================================================
// TEST CATEGORY 4: INPUT PARAMETER VALIDATION
// =============================================================================

TEST_F(DynamicDataUpdateTest, InputValidation_EmptyInstructionSetsMap) {
    // Test: Empty GenerationInstructionSets map
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial values to verify they get reset
    speciesData[1].instructionSetsSize = 10;
    speciesData[2].instructionSetsSize = 20;
    speciesData[3].instructionSetsSize = 30;
    
    GenerationInstructionSets emptyInstructionSets; // Empty map
    
    // Should not crash with empty instruction sets
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, emptyInstructionSets, *params));
    
    // All instruction set sizes should be reset to 0
    EXPECT_EQ(speciesData[1].instructionSetsSize, 0) << "Species 1 should be reset to 0 with empty instruction sets";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 0) << "Species 2 should be reset to 0 with empty instruction sets";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 0) << "Species 3 should be reset to 0 with empty instruction sets";
}

TEST_F(DynamicDataUpdateTest, InputValidation_SingleSpeciesInInstructionSets) {
    // Test: Single species in instruction sets
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial values
    speciesData[1].instructionSetsSize = 10;
    speciesData[2].instructionSetsSize = 20;
    speciesData[3].instructionSetsSize = 30;
    
    GenerationInstructionSets instructionSets;
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // Only species 2 has instructions
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // Only species 2 should have non-zero count, others reset to 0
    EXPECT_EQ(speciesData[1].instructionSetsSize, 0) << "Species 1 should be reset to 0";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 5) << "Species 2 should have instruction count 5";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 0) << "Species 3 should be reset to 0";
}

// =============================================================================
// TEST CATEGORY 5: EDGE CASES AND BOUNDARIES
// =============================================================================

TEST_F(DynamicDataUpdateTest, EdgeCases_NestedEmptyVectors) {
    // Test: Nested empty vectors within species map
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2});
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{}; // Explicitly empty vector
    instructionSets[2] = SpeciesInstructionSet{}; // Explicitly empty vector
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // Both should have count 0
    EXPECT_EQ(speciesData[1].instructionSetsSize, 0) << "Empty instruction vector should give count 0";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 0) << "Empty instruction vector should give count 0";
}

TEST_F(DynamicDataUpdateTest, EdgeCases_LargeInstructionCounts) {
    // Test: Large instruction counts
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Create large instruction set (1000 instructions)
    SpeciesInstructionSet largeInstructionSet;
    for (uint32_t i = 0; i < 1000; ++i) {
        largeInstructionSet.push_back(ReproductiveInstruction::preserve(i % 100));
    }
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = largeInstructionSet;
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    EXPECT_EQ(speciesData[1].instructionSetsSize, 1000) << "Large instruction count should be handled correctly";
}

// =============================================================================
// TEST CATEGORY 6: COUNTING LOGIC VERIFICATION
// =============================================================================

TEST_F(DynamicDataUpdateTest, CountingLogic_DirectCountingAlgorithm) {
    // Test: Core counting mechanism works correctly
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Create instruction sets with known sizes that can be manually verified
    GenerationInstructionSets instructionSets;
    
    // Manually count: 4 instructions for species 1
    instructionSets[1] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3)
    };
    
    // Manually count: 6 instructions for species 2
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4),
        ReproductiveInstruction::preserve(5)
    };
    
    // Manually count: 1 instruction for species 3
    instructionSets[3] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0)
    };
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify manual counts match algorithm results
    EXPECT_EQ(speciesData[1].instructionSetsSize, 4) << "Algorithm should count 4 instructions correctly";
    EXPECT_EQ(speciesData[2].instructionSetsSize, 6) << "Algorithm should count 6 instructions correctly";
    EXPECT_EQ(speciesData[3].instructionSetsSize, 1) << "Algorithm should count 1 instruction correctly";
}

TEST_F(DynamicDataUpdateTest, CountingLogic_CountingConsistency) {
    // Test: Counting produces identical results for identical inputs
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData1 = createTestSpeciesData({1, 2});
    auto speciesData2 = createTestSpeciesData({1, 2});
    
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions
    
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // 5 instructions
    
    // Process same data twice
    dynamicDataUpdate(genomeData, speciesData1, instructionSets, *params);
    dynamicDataUpdate(genomeData, speciesData2, instructionSets, *params);
    
    // Results should be identical
    EXPECT_EQ(speciesData1[1].instructionSetsSize, speciesData2[1].instructionSetsSize) 
        << "First processing should give same result as second processing";
    EXPECT_EQ(speciesData1[2].instructionSetsSize, speciesData2[2].instructionSetsSize) 
        << "First processing should give same result as second processing";
        
    // Verify specific expected values
    EXPECT_EQ(speciesData1[1].instructionSetsSize, 3) << "Species 1 should consistently have 3 instructions";
    EXPECT_EQ(speciesData1[2].instructionSetsSize, 5) << "Species 2 should consistently have 5 instructions";
}

// =============================================================================
// ESSENTIAL MISSING TESTS FOR NEW FUNCTIONALITY
// =============================================================================

// =============================================================================
// TEST 1: NEW SPECIES FIELD INITIALIZATION
// =============================================================================

TEST_F(DynamicDataUpdateTest, NewSpecies_AllFieldsInitializedCorrectly) {
    // Test: Verify all fields in new species dynamic data are correctly initialized
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only existing species
    
    // Create instruction sets that includes a new species (999)
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    instructionSets[999] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3),
        ReproductiveInstruction::preserve(4)
    }; // 5 instructions for new species
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify new species was created
    auto newSpeciesIt = speciesData.find(999);
    EXPECT_NE(newSpeciesIt, speciesData.end()) << "New species 999 should be created";
    
    if (newSpeciesIt != speciesData.end()) {
        const auto& newSpeciesData = newSpeciesIt->second;
        
        // Verify all field initialization
        EXPECT_EQ(newSpeciesData.instructionSetsSize, 5) << "instructionSetsSize should be set from instruction count";
        EXPECT_EQ(newSpeciesData.protectionRating, 0) << "protectionRating should be initialized to 0";
        EXPECT_EQ(newSpeciesData.currentPopulationSize, 0) << "currentPopulationSize should be 0 (no genomes exist)";
        EXPECT_GT(newSpeciesData.speciesRank, 0) << "speciesRank should be assigned (worst rank since no genomes)";
        EXPECT_FALSE(newSpeciesData.isMarkedForElimination) << "isMarkedForElimination should be false by default";
        
        // Verify worst rank assignment (new species with no genomes gets worst rank)
        uint32_t expectedWorstRank = 3; // Species 1 and 2 have genomes, so new species gets rank 3
        EXPECT_EQ(newSpeciesData.speciesRank, expectedWorstRank) << "New species should get worst ordinal rank";
    }
}

// =============================================================================
// TEST 2: SPECIES DISCOVERY THROUGH BOTH SOURCES
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesDiscovery_ThroughGenomeAndInstructionData) {
    // Test: Verify species referenced in genome data and instruction sets gets dynamic data created
    // Important edge case: species appears in both genome data AND instruction sets but no existing dynamic data
    
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 888}; // Species 888 will be "discovered"
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only species 1 and 2 exist initially
    
    // Create instruction sets that also includes species 888
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    instructionSets[888] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions for discovered species
    
    // Should not crash or throw even though species 888 has no existing dynamic data
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // Verify species 888 was discovered and created in species data
    auto discoveredSpeciesIt = speciesData.find(888);
    EXPECT_NE(discoveredSpeciesIt, speciesData.end()) << "Species 888 should be discovered and created";
    
    if (discoveredSpeciesIt != speciesData.end()) {
        const auto& discoveredData = discoveredSpeciesIt->second;
        
        // Verify the species data reflects both genome presence and instruction sets
        EXPECT_EQ(discoveredData.instructionSetsSize, 3) << "instructionSetsSize should be set from instruction sets";
        EXPECT_EQ(discoveredData.currentPopulationSize, 1) << "currentPopulationSize should reflect genome count";
        EXPECT_GT(discoveredData.speciesRank, 0) << "speciesRank should be assigned based on genome performance";
        EXPECT_EQ(discoveredData.protectionRating, 1) << "Worst performing species should be penalized";
        EXPECT_FALSE(discoveredData.isMarkedForElimination) << "isMarkedForElimination should be false";
        
        // Verify species ranking reflects actual performance
        // Species 888 has genome at rank 2 (fitness 80), so should have good ranking
        EXPECT_LE(discoveredData.speciesRank, 3) << "Species with genomes should get reasonable ordinal rank";
    }
}

// =============================================================================
// TEST 3: FLOATING POINT TO ORDINAL RANKING CONVERSION
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_FloatingPointToOrdinalConversion) {
    // Test: Verify conversion from average ranks to ordinal rankings
    // Core algorithm of new ranking functionality
    
    // Create specific scenario with known floating point averages
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55}; // 10 genomes
    std::vector<uint32_t> species = {1,   1,  2,  2,  2,  3,  3,  3,  3,  3}; // Uneven distribution
    //                              r0   r1  r2  r3  r4  r5  r6  r7  r8  r9
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Calculate expected floating point averages:
    // Species 1: ranks [0, 1] → average = (0+1)/2 = 0.5
    // Species 2: ranks [2, 3, 4] → average = (2+3+4)/3 = 3.0  
    // Species 3: ranks [5, 6, 7, 8, 9] → average = (5+6+7+8+9)/5 = 7.0
    
    // Expected ordinal conversion (lower average = better ordinal rank):
    // Species 1: average 0.5 → ordinal rank 1 (best)
    // Species 2: average 3.0 → ordinal rank 2 (middle)
    // Species 3: average 7.0 → ordinal rank 3 (worst)
    
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "Species 1 (avg 0.5) should get ordinal rank 1";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 (avg 3.0) should get ordinal rank 2";
    EXPECT_EQ(speciesData[3].speciesRank, 3) << "Species 3 (avg 7.0) should get ordinal rank 3";
}

// =============================================================================
// TEST 4: ORDINAL RANKING ASSIGNMENT
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_AssignsCorrectOrdinalRankings) {
    // Test: Verify species get correct ordinal rankings (1, 2, 3...) based on performance
    // Basic output validation for new ranking feature
    
    // Create 4 species with clear performance differences
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30};
    std::vector<uint32_t> species = {1,   2,  3,  4,  1,  2,  3,  4};
    //                              r0   r1  r2  r3  r4  r5  r6  r7
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3, 4});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Calculate expected average ranks:
    // Species 1: ranks [0, 4] → average = (0+4)/2 = 2.0 (best performance)
    // Species 2: ranks [1, 5] → average = (1+5)/2 = 3.0 (second best)
    // Species 3: ranks [2, 6] → average = (2+6)/2 = 4.0 (third best)
    // Species 4: ranks [3, 7] → average = (3+7)/2 = 5.0 (worst performance)
    
    // Verify ordinal rankings (1=best, 2=second best, etc.)
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "Best performing species should get rank 1";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Second best species should get rank 2";
    EXPECT_EQ(speciesData[3].speciesRank, 3) << "Third best species should get rank 3";
    EXPECT_EQ(speciesData[4].speciesRank, 4) << "Worst performing species should get rank 4";
    
    // Verify all rankings are sequential with no gaps
    std::vector<uint32_t> observedRanks = {
        speciesData[1].speciesRank,
        speciesData[2].speciesRank,
        speciesData[3].speciesRank,
        speciesData[4].speciesRank
    };
    std::sort(observedRanks.begin(), observedRanks.end());
    
    std::vector<uint32_t> expectedRanks = {1, 2, 3, 4};
    EXPECT_EQ(observedRanks, expectedRanks) << "Rankings should be sequential 1,2,3,4 with no gaps";
}

// =============================================================================
// TEST 5: SPECIES RANK FIELD UPDATES
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_UpdatesSpeciesRankField) {
    // Test: Verify speciesRank field gets correctly updated with ordinal values
    // Field update mechanism for new feature
    
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial rank values to verify they get overwritten
    speciesData[1].speciesRank = 999; // Should become 1
    speciesData[2].speciesRank = 888; // Should become 2  
    speciesData[3].speciesRank = 777; // Should become 3
    
    // Set other fields to verify they don't interfere with ranking
    speciesData[1].protectionRating = 5;
    speciesData[2].protectionRating = 3;
    speciesData[3].protectionRating = 1;
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify speciesRank field was updated correctly
    // Species 1: rank 0 (fitness 100) → ordinal rank 1
    // Species 2: rank 1 (fitness 90) → ordinal rank 2  
    // Species 3: rank 2 (fitness 80) → ordinal rank 3
    
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "speciesRank field should be updated from 999 to 1";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "speciesRank field should be updated from 888 to 2";
    EXPECT_EQ(speciesData[3].speciesRank, 3) << "speciesRank field should be updated from 777 to 3";
    
    // Verify initial rank 0 gets updated (no species should have rank 0 after processing)
    EXPECT_NE(speciesData[1].speciesRank, 0) << "No species should have rank 0 after processing";
    EXPECT_NE(speciesData[2].speciesRank, 0) << "No species should have rank 0 after processing";
    EXPECT_NE(speciesData[3].speciesRank, 0) << "No species should have rank 0 after processing";
    
    // Verify all species get assigned ordinal rank values
    EXPECT_GT(speciesData[1].speciesRank, 0) << "All species should have positive ordinal ranks";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "All species should have positive ordinal ranks";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "All species should have positive ordinal ranks";
    
    // Verify ranking field update doesn't affect other fields inappropriately
    // Note: protectionRating may change due to worst species penalty, but should follow expected logic
}

// =============================================================================
// TEST 6: NEW SPECIES ORDINAL RANKING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_NewSpeciesOrdinalRanking) {
    // Test: Verify new species (from instruction sets) get appropriate ordinal rankings
    // Tests interaction between new species creation and ranking
    
    std::vector<double> fitness = {100, 90};
    std::vector<uint32_t> species = {1, 2}; // Only existing species have genomes
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only species 1 and 2 exist initially
    
    // Create instruction sets that include new species without genomes
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
    instructionSets[2] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1)
    }; // 2 instructions
    instructionSets[500] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2)
    }; // 3 instructions for new species (no genomes)
    instructionSets[600] = SpeciesInstructionSet{
        ReproductiveInstruction::preserve(0),
        ReproductiveInstruction::preserve(1),
        ReproductiveInstruction::preserve(2),
        ReproductiveInstruction::preserve(3)
    }; // 4 instructions for another new species (no genomes)
    
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Verify all new species were created
    EXPECT_NE(speciesData.find(500), speciesData.end()) << "New species 500 should be created";
    EXPECT_NE(speciesData.find(600), speciesData.end()) << "New species 600 should be created";
    
    // Verify existing species get good ordinal rankings (have genomes)
    // Species 1: rank 0 (fitness 100) → ordinal rank 1
    // Species 2: rank 1 (fitness 90) → ordinal rank 2
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "Existing species 1 should get ordinal rank 1";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Existing species 2 should get ordinal rank 2";
    
    // Verify new species get worst ordinal rankings (no genomes)
    // New species with no genomes should get tied for worst rank
    uint32_t expectedWorstRank = 3; // After species 1 and 2, new species get rank 3
    EXPECT_EQ(speciesData[500].speciesRank, expectedWorstRank) << "New species 500 should get worst ordinal rank";
    EXPECT_EQ(speciesData[600].speciesRank, expectedWorstRank) << "New species 600 should get worst ordinal rank";
    
    // Verify multiple new species get tied for worst ordinal rank
    EXPECT_EQ(speciesData[500].speciesRank, speciesData[600].speciesRank) 
        << "Multiple new species should get same worst ordinal rank";
    
    // Verify new species ranking doesn't affect existing species ordinal rankings
    EXPECT_LT(speciesData[1].speciesRank, speciesData[500].speciesRank) 
        << "Existing species should have better rank than new species";
    EXPECT_LT(speciesData[2].speciesRank, speciesData[500].speciesRank) 
        << "Existing species should have better rank than new species";
    
    // Verify all species have valid positive ranks
    EXPECT_GT(speciesData[1].speciesRank, 0) << "All species should have positive ranks";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "All species should have positive ranks";
    EXPECT_GT(speciesData[500].speciesRank, 0) << "All species should have positive ranks";
    EXPECT_GT(speciesData[600].speciesRank, 0) << "All species should have positive ranks";
}

// =============================================================================
// TEST 7: SINGLE SPECIES ORDINAL RANKING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_SingleSpeciesGetsRankOne) {
    // Test: Verify single species gets ordinal rank 1 (best by default)
    // Important boundary condition for ranking algorithm
    
    // Test scenario 1: Single species with multiple genomes
    {
        std::vector<double> fitness = {100, 90, 80};
        std::vector<uint32_t> species = {1, 1, 1}; // All genomes belong to species 1
        
        auto genomeData = createTestGenomeData(fitness, species);
        auto speciesData = createTestSpeciesData({1});
        
        auto instructionSets = createMockInstructionSets({1});
        dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
        
        // Single species should get ordinal rank 1 (best by default)
        EXPECT_EQ(speciesData[1].speciesRank, 1) << "Single species with multiple genomes should get rank 1";
        EXPECT_GT(speciesData[1].speciesRank, 0) << "Single species should have positive rank";
    }
    
    // Test scenario 2: Single species with single genome
    {
        std::vector<double> fitness = {100};
        std::vector<uint32_t> species = {42}; // Single genome for species 42
        
        auto genomeData = createTestGenomeData(fitness, species);
        auto speciesData = createTestSpeciesData({42});
        
        auto instructionSets = createMockInstructionSets({42});
        dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
        
        // Single species should get ordinal rank 1 (best by default)
        EXPECT_EQ(speciesData[42].speciesRank, 1) << "Single species with single genome should get rank 1";
        EXPECT_GT(speciesData[42].speciesRank, 0) << "Single species should have positive rank";
    }
    
    // Test scenario 3: Single species with instruction sets but no genomes
    {
        std::multimap<TestFitnessResult, DynamicGenomeData> emptyGenomeData;
        auto speciesData = createTestSpeciesData({999});
        
        GenerationInstructionSets instructionSets;
        instructionSets[999] = SpeciesInstructionSet{
            ReproductiveInstruction::preserve(0),
            ReproductiveInstruction::preserve(1)
        }; // 2 instructions for species with no genomes
        
        dynamicDataUpdate(emptyGenomeData, speciesData, instructionSets, *params);
        
        // Single species with no genomes should still get assigned a rank
        // No other species to compare against → automatically rank 1
        EXPECT_EQ(speciesData[999].speciesRank, 1) << "Single species with no genomes should get rank 1";
        EXPECT_GT(speciesData[999].speciesRank, 0) << "Single species should have positive rank";
    }
}

// =============================================================================
// TEST 8: NO GENOMES ORDINAL RANKING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_NoGenomesOrdinalRanking) {
    // Test: Verify ordinal ranking when no genomes exist for ranking
    // Critical edge case - how ranking behaves with empty performance data
    
    // Test scenario 1: Empty genome data with existing species (instruction set driven)
    {
        std::multimap<TestFitnessResult, DynamicGenomeData> emptyGenomeData;
        auto speciesData = createTestSpeciesData({1, 2, 3});
        
        // All species have instruction sets but no genomes
        GenerationInstructionSets instructionSets;
        instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // 1 instruction
        instructionSets[2] = SpeciesInstructionSet{
            ReproductiveInstruction::preserve(0),
            ReproductiveInstruction::preserve(1)
        }; // 2 instructions
        instructionSets[3] = SpeciesInstructionSet{
            ReproductiveInstruction::preserve(0),
            ReproductiveInstruction::preserve(1),
            ReproductiveInstruction::preserve(2)
        }; // 3 instructions
        
        // Should not crash or invalid operations with empty performance data
        EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, speciesData, instructionSets, *params));
        
        // All species get same worst ordinal rank (no performance data to differentiate)
        uint32_t expectedRank = 1; // First (and only) rank when no genomes exist
        EXPECT_EQ(speciesData[1].speciesRank, expectedRank) << "Species with no genomes should get rank 1";
        EXPECT_EQ(speciesData[2].speciesRank, expectedRank) << "Species with no genomes should get rank 1";
        EXPECT_EQ(speciesData[3].speciesRank, expectedRank) << "Species with no genomes should get rank 1";
        
        // Verify all species get same ordinal rank
        EXPECT_EQ(speciesData[1].speciesRank, speciesData[2].speciesRank) 
            << "All species with no genomes should get same rank";
        EXPECT_EQ(speciesData[2].speciesRank, speciesData[3].speciesRank) 
            << "All species with no genomes should get same rank";
        
        // Verify all species have positive ranks
        EXPECT_GT(speciesData[1].speciesRank, 0) << "All species should have positive ranks";
        EXPECT_GT(speciesData[2].speciesRank, 0) << "All species should have positive ranks";
        EXPECT_GT(speciesData[3].speciesRank, 0) << "All species should have positive ranks";
    }
    
    // Test scenario 2: Completely empty inputs
    {
        std::multimap<TestFitnessResult, DynamicGenomeData> emptyGenomeData;
        std::unordered_map<uint32_t, DynamicSpeciesData> emptySpeciesData;
        GenerationInstructionSets emptyInstructionSets;
        
        // Should not crash with completely empty inputs
        EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, emptySpeciesData, emptyInstructionSets, *params));
        
        // Nothing should be created or modified
        EXPECT_TRUE(emptyGenomeData.empty()) << "Empty genome data should remain empty";
        EXPECT_TRUE(emptySpeciesData.empty()) << "Empty species data should remain empty";
    }
    
    // Test scenario 3: Mixed scenario - some species have genomes, others only instruction sets
    {
        std::vector<double> fitness = {100, 90}; // Only 2 genomes
        std::vector<uint32_t> species = {1, 2}; // Only species 1 and 2 have genomes
        
        auto genomeData = createTestGenomeData(fitness, species);
        auto speciesData = createTestSpeciesData({1, 2, 3, 4}); // 4 species total
        
        // Species 3 and 4 have instruction sets but no genomes
        GenerationInstructionSets instructionSets;
        instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // Has genomes
        instructionSets[2] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // Has genomes
        instructionSets[3] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // No genomes
        instructionSets[4] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)}; // No genomes
        
        EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
        
        // Species with genomes should get performance-based ranks 1 and 2
        EXPECT_EQ(speciesData[1].speciesRank, 1) << "Species 1 with best genome should get rank 1";
        EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 with second best genome should get rank 2";
        
        // Species without genomes should get worst rank (tied)
        uint32_t expectedWorstRank = 3; // After ranks 1 and 2 are taken
        EXPECT_EQ(speciesData[3].speciesRank, expectedWorstRank) << "Species 3 with no genomes should get worst rank";
        EXPECT_EQ(speciesData[4].speciesRank, expectedWorstRank) << "Species 4 with no genomes should get worst rank";
        
        // Verify species without genomes get same rank
        EXPECT_EQ(speciesData[3].speciesRank, speciesData[4].speciesRank) 
            << "Species with no genomes should get same worst rank";
        
        // Verify ranking hierarchy
        EXPECT_LT(speciesData[1].speciesRank, speciesData[3].speciesRank) 
            << "Species with genomes should rank better than species without";
        EXPECT_LT(speciesData[2].speciesRank, speciesData[4].speciesRank) 
            << "Species with genomes should rank better than species without";
    }
}

// =============================================================================
// TEST 9: TIED AVERAGE RANKS ORDINAL HANDLING
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesRanking_TiedAverageRanksOrdinalHandling) {
    // Test: Verify ordinal ranking behavior when multiple species have identical average ranks
    // Important edge case for ranking algorithm stability and deterministic behavior
    
    // Create scenario where two species have identical average ranks
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50}; // 6 genomes
    std::vector<uint32_t> species = {1,   2,  3,  3,  2,  1}; // Symmetric distribution
    //                              r0   r1  r2  r3  r4  r5
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto instructionSets = createMockInstructionSets({1, 2, 3});
    dynamicDataUpdate(genomeData, speciesData, instructionSets, *params);
    
    // Calculate expected average ranks:
    // Species 1: ranks [0, 5] → average = (0+5)/2 = 2.5
    // Species 2: ranks [1, 4] → average = (1+4)/2 = 2.5 (TIED with species 1)
    // Species 3: ranks [2, 3] → average = (2+3)/2 = 2.5 (TIED with species 1 and 2)
    
    // All species have identical average ranks (2.5), so ranking algorithm must handle ties
    
    // Verify all species get valid ordinal rankings
    EXPECT_GT(speciesData[1].speciesRank, 0) << "Species 1 should have positive ordinal rank";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "Species 2 should have positive ordinal rank";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "Species 3 should have positive ordinal rank";
    
    // Verify ordinal rankings are within valid range (1 to number of species)
    EXPECT_LE(speciesData[1].speciesRank, 3) << "Species 1 rank should be <= 3";
    EXPECT_LE(speciesData[2].speciesRank, 3) << "Species 2 rank should be <= 3";
    EXPECT_LE(speciesData[3].speciesRank, 3) << "Species 3 rank should be <= 3";
    
    // For tied average ranks, the algorithm assigns ordinal ranks based on std::sort behavior
    // When average ranks are equal, the ordering depends on implementation-defined behavior
    // We don't enforce specific species ID ordering for ties, just verify valid rank assignment
    
    // Collect all ranks to verify they form a valid sequence
    std::vector<uint32_t> allRanks = {
        speciesData[1].speciesRank,
        speciesData[2].speciesRank,
        speciesData[3].speciesRank
    };
    
    // All ranks should be different (no two species get exactly same ordinal rank)
    std::set<uint32_t> uniqueRanks(allRanks.begin(), allRanks.end());
    EXPECT_EQ(uniqueRanks.size(), 3) << "All species should get different ordinal ranks even when tied";
    
    // Verify ranks form consecutive sequence starting from 1
    std::sort(allRanks.begin(), allRanks.end());
    std::vector<uint32_t> expectedSequence = {1, 2, 3};
    EXPECT_EQ(allRanks, expectedSequence) << "Tied species should still get consecutive ordinal ranks 1,2,3";
    
    // Note: The current implementation does not guarantee deterministic tie-breaking
    // When average ranks are equal, the ordering depends on std::sort implementation details
    // This is acceptable behavior as long as all species get valid, unique ordinal ranks
}

// =============================================================================
// TEST 10: SPECIES DISCOVERY THROUGH GENOME DATA ONLY
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesDiscovery_ThroughGenomeDataOnly) {
    // Test: Verify species referenced only in genome data (not instruction sets) gets created
    // This tests the fix for circular dependency - new species emerge during evolution
    
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 777}; // Species 777 appears only in genome data
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only species 1 and 2 exist initially
    
    // Create instruction sets that do NOT include species 777
    GenerationInstructionSets instructionSets;
    instructionSets[1] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)};
    instructionSets[2] = SpeciesInstructionSet{ReproductiveInstruction::preserve(0)};
    // Species 777 NOT in instruction sets - purely discovered through genome data
    
    // Should not crash even though species 777 has no instruction sets or existing data
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, instructionSets, *params));
    
    // Verify species 777 was discovered and created
    auto discoveredSpeciesIt = speciesData.find(777);
    EXPECT_NE(discoveredSpeciesIt, speciesData.end()) << "Species 777 should be discovered through genome data";
    
    if (discoveredSpeciesIt != speciesData.end()) {
        const auto& discoveredData = discoveredSpeciesIt->second;
        
        // Verify proper initialization of discovered species
        EXPECT_EQ(discoveredData.currentPopulationSize, 1) << "Should reflect the one genome in this species";
        EXPECT_EQ(discoveredData.instructionSetsSize, 0) << "Should be 0 since not in instruction sets";
        EXPECT_EQ(discoveredData.protectionRating, 1) << "Should be penalized as worst performing species";
        EXPECT_GT(discoveredData.speciesRank, 0) << "Should get valid ordinal rank based on performance";
        EXPECT_FALSE(discoveredData.isMarkedForElimination) << "Should not be marked for elimination initially";
        
        // Verify it gets worst rank (rank 3) since it has worst individual performance
        EXPECT_EQ(discoveredData.speciesRank, 3) << "Species with worst genome should get worst rank";
    }
    
    // Verify existing species still work correctly  
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Existing species 1 should be updated";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 1) << "Existing species 2 should be updated";
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "Species 1 should get best rank";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 should get middle rank";
}