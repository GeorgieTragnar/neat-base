#include <gtest/gtest.h>
#include <map>
#include <unordered_map>

#include "tests/test_common.h"
#include "version3/population/DynamicDataUpdate.hpp"

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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, zeroParams);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, hundredParams);
    
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
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, multiWorstParams);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params); // worstSpeciesCount = 1
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params); // worstSpeciesCount = 1
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params); // worstSpeciesCount = 1
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, multiWorstParams);
    
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
    
    dynamicDataUpdate(genomeData1, speciesData, *params);
    
    // Species 3 should increment to 2
    EXPECT_EQ(speciesData[3].protectionRating, 2) << "Species 3 should increment to 2 in first generation";
    
    // Second generation: Species 3 improves, Species 2 becomes worst
    std::vector<double> fitness2 = {100, 95, 50, 45, 90, 85}; // Species 3 now has better fitness
    std::vector<uint32_t> species2 = {1,   1,  2,  2,  3,  3};
    
    auto genomeData2 = createTestGenomeData(fitness2, species2);
    
    dynamicDataUpdate(genomeData2, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, emptySpeciesData, *params));
    
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
    
    EXPECT_NO_THROW(dynamicDataUpdate(emptyGenomeData, speciesData, *params));
    
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
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, emptySpeciesData, *params));
    
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
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, excessParams));
    
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
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, equalParams));
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    dynamicDataUpdate(genomeData, speciesData, *params);
    
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
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, *params));
    
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
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, highProtectionParams));
    
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
    
    dynamicDataUpdate(genomeData, speciesData, fullProtectionParams);
    
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
    
    EXPECT_NO_THROW(dynamicDataUpdate(genomeData, speciesData, *params));
    
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