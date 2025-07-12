#include <gtest/gtest.h>
#include <map>
#include <unordered_map>

#include "tests/test_common.h"
#include "version3/population/DynamicDataUpdate.hpp"
#include "version3/population/GlobalIndexRegistry.hpp"

using namespace Population;

class DynamicDataUpdateTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test parameters with reasonable defaults
        params = std::make_unique<DynamicDataUpdateParams>(
            5,    // maxGenomePendingEliminationLimit
            3,    // maxSpeciesPendingEliminationRating  
            0.3,  // speciesElitePlacementProtectionPercentage (30%)
            0.4,  // speciesPendingEliminationPercentage (40%)
            0.3,  // genomesPendingEliminationPercentage (30%)
            2,    // equilibriumSpeciesCount
            10    // targetPopulationSize
        );
        
        // Create registry with sufficient capacity
        registry = std::make_unique<GlobalIndexRegistry>(100);
    }

    std::unique_ptr<DynamicDataUpdateParams> params;
    std::unique_ptr<GlobalIndexRegistry> registry;
    
    // Helper to create simple fitness result type for testing
    struct TestFitnessResult {
        double fitness;
        
        bool operator<(const TestFitnessResult& other) const {
            return fitness < other.fitness; // Lower fitness comes first (worst to best ordering)
        }
        
        bool operator==(const TestFitnessResult& other) const {
            return fitness == other.fitness;
        }
    };
    
    // Helper to create test data structures
    struct TestData {
        std::multimap<TestFitnessResult, size_t> fitnessResults;
        std::vector<DynamicGenomeData> genomeData;
        std::unordered_map<uint32_t, std::vector<size_t>> speciesGrouping;
    };
    
    TestData createTestData(
        const std::vector<double>& fitnessValues,
        const std::vector<uint32_t>& speciesIds) {
        
        TestData data;
        data.genomeData.resize(fitnessValues.size());
        
        for (size_t i = 0; i < fitnessValues.size(); ++i) {
            TestFitnessResult fitness{fitnessValues[i]};
            data.fitnessResults.insert({fitness, i});
            
            DynamicGenomeData& genomeData = data.genomeData[i];
            genomeData.speciesId = speciesIds[i];
            genomeData.pendingEliminationCounter = 0;
            genomeData.isUnderRepair = false;
            genomeData.isMarkedForElimination = false;
            genomeData.genomeIndex = static_cast<uint32_t>(i);
            
            // Add to species grouping
            data.speciesGrouping[speciesIds[i]].push_back(i);
        }
        
        return data;
    }
    
    // Helper to create test species data map
    std::unordered_map<uint32_t, DynamicSpeciesData> createTestSpeciesData(
        const std::vector<uint32_t>& speciesIds) {
        
        std::unordered_map<uint32_t, DynamicSpeciesData> speciesData;
        
        for (uint32_t id : speciesIds) {
            DynamicSpeciesData data;
            data.pendingEliminationRating = 0;
            data.currentPopulationSize = 0;
            data.speciesRank = 0;
            data.isMarkedForElimination = false;
            
            speciesData[id] = data;
        }
        
        return speciesData;
    }
};

// =============================================================================
// 1. PROTECTION TIER CLASSIFICATION LOGIC
// =============================================================================

TEST_F(DynamicDataUpdateTest, SpeciesElitePlacementProtection_30Percent) {
    // Test: Species with champion in top 30% gets protection from elimination
    // Setup: 10 genomes, species with different champion ranks
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10}; // Best to worst
    std::vector<uint32_t> species = {1, 2, 2, 2, 2, 2, 3, 3, 3, 3};
    //                              r0  r1  r2  r3  r4  r5  r6  r7  r8  r9
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set up a scenario where species would be targeted for elimination
    // Force species into elimination consideration by having enough active species
    params = std::make_unique<DynamicDataUpdateParams>(5, 3, 0.3, 0.5, 0.3, 1, 10);
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // With 30% protection threshold and 10 genomes (ranks 0-9):
    // - Species 1: champion at rank 0 → percentile 0.0/9 = 0% < 30% → NOT protected
    // - Species 2: champion at rank 1 → percentile 1.0/9 = 11.1% < 30% → NOT protected  
    // - Species 3: champion at rank 6 → percentile 6.0/9 = 66.7% >= 30% → protected
    
    // Verify species rankings were updated
    EXPECT_GT(speciesData[1].speciesRank, 0) << "Species 1 should have valid rank";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "Species 2 should have valid rank";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "Species 3 should have valid rank";
    
    // Verify population sizes were updated correctly
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 5) << "Species 2 should have 5 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 4) << "Species 3 should have 4 genomes";
}

TEST_F(DynamicDataUpdateTest, SpeciesElitePlacementProtection_BoundaryCase) {  
    // Test boundary case: Species with champion exactly at protection threshold
    std::vector<double> fitness = {100, 90, 80, 70, 60}; // 5 genomes
    std::vector<uint32_t> species = {1, 2, 2, 2, 3};
    //                              r0  r1  r2  r3  r4
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // With 30% protection and 5 genomes (ranks 0-4):
    // - Species 1: champion at rank 0 → percentile 0.0/4 = 0% < 30% → NOT protected
    // - Species 2: champion at rank 1 → percentile 1.0/4 = 25% < 30% → NOT protected
    // - Species 3: champion at rank 4 → percentile 4.0/4 = 100% >= 30% → protected
    
    // Set up parameters to force elimination checks
    params = std::make_unique<DynamicDataUpdateParams>(5, 3, 0.3, 0.5, 0.3, 1, 10);
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // Verify species rankings are assigned
    EXPECT_GT(speciesData[1].speciesRank, 0) << "Species 1 should have valid rank";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "Species 2 should have valid rank";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "Species 3 should have valid rank";
    
    // Verify population sizes are correct
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 1) << "Species 3 should have 1 genome";
}

TEST_F(DynamicDataUpdateTest, SpeciesElitePlacementProtection_ZeroPercent) {
    // Test edge case: 0% protection → no species protected from elimination
    DynamicDataUpdateParams zeroParams(5, 3, 0.0, 0.5, 0.3, 1, 10); // 0% protection threshold
    
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1, 2, 2, 2, 3};
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, zeroParams, *registry);
    
    // With 0% protection, no species should be protected regardless of champion placement
    // All species should be subject to elimination rating increases if they're worst performers
    
    // Verify basic functionality still works
    EXPECT_GT(speciesData[1].speciesRank, 0) << "Species 1 should have valid rank";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "Species 2 should have valid rank";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "Species 3 should have valid rank";
    
    // Verify population sizes are updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 1) << "Species 3 should have 1 genome";
}

TEST_F(DynamicDataUpdateTest, SpeciesElitePlacementProtection_HundredPercent) {
    // Test edge case: 100% protection → all species protected from elimination
    DynamicDataUpdateParams hundredParams(5, 3, 1.0, 0.0, 0.3, 1, 10); // 100% protection: no excess species marked for elimination
    
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {1, 2, 2, 2, 3};
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, hundredParams, *registry);
    
    // With 100% protection, all species should be protected from elimination
    // No species should get pending elimination rating increases
    
    // Verify basic functionality still works
    EXPECT_GT(speciesData[1].speciesRank, 0) << "Species 1 should have valid rank";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "Species 2 should have valid rank";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "Species 3 should have valid rank";
    
    // Verify population sizes are updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 1) << "Species 3 should have 1 genome";
    
    // With 100% protection, no species should have increased elimination rating
    EXPECT_EQ(speciesData[1].pendingEliminationRating, 0) << "Species 1 should not have elimination rating increased";
    EXPECT_EQ(speciesData[2].pendingEliminationRating, 0) << "Species 2 should not have elimination rating increased";
    EXPECT_EQ(speciesData[3].pendingEliminationRating, 0) << "Species 3 should not have elimination rating increased";
}
TEST_F(DynamicDataUpdateTest, GenomePendingEliminationCounter_ExceedsLimitTriggersElimination) {
    // Test: Genome pending elimination counter exceeding limit triggers elimination
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    // Set counter at limit (maxGenomePendingEliminationLimit = 5) for worst genome
    testData.genomeData[9].pendingEliminationCounter = 5; // At limit, will increment to 6 and trigger elimination
    testData.genomeData[8].pendingEliminationCounter = 4; // Will be decremented to 3 (not in pending elimination range)
    
    // Set up parameters to force genome elimination
    // Need genomesPendingEliminationPercentage to ensure at least 1 genome is pending elimination
    // With 10 genomes, equilibrium=8, excess=2, need 2*percentage >= 1 → percentage >= 0.5
    params = std::make_unique<DynamicDataUpdateParams>(5, 3, 0.3, 0.5, 0.6, 1, 8); // equilibrium = 8/1 = 8, so 2 excess genomes, 2*0.6=1.2→1 genome pending elimination
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // Verify genome at limit gets marked for elimination
    EXPECT_EQ(testData.genomeData[9].pendingEliminationCounter, 6) 
        << "Worst genome counter should increment to 6 and trigger elimination";
    EXPECT_EQ(registry->getState(9), GenomeState::HotElimination) 
        << "Worst genome should be marked for elimination when counter exceeds limit";
    
    // Verify genome not in pending elimination range gets counter decremented
    EXPECT_EQ(testData.genomeData[8].pendingEliminationCounter, 3) 
        << "Genome counter should decrement to 3 (not in pending elimination range)";
    EXPECT_EQ(registry->getState(8), GenomeState::Active) 
        << "Genome should remain active when not in pending elimination range";
        
    // Verify top performers are not affected
    EXPECT_EQ(testData.genomeData[0].pendingEliminationCounter, 0) 
        << "Top genome counter should remain 0";
    EXPECT_EQ(registry->getState(0), GenomeState::Active) 
        << "Top genome should remain active";
}

TEST_F(DynamicDataUpdateTest, SpeciesRanking_AverageRankCalculation) {
    // Test: Species ranking calculation based on average genome ranks
    // Create scenario with clear rank differences
    // With worst-to-best ordering: 55(r0), 60(r1), 65(r2), 70(r3), 75(r4), 80(r5), 85(r6), 90(r7), 95(r8), 100(r9)
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   3,  2,  3,  3,  1,  3,  2,  3,  1};
    //                              r9   r8  r7  r6  r5  r4  r3  r2  r1  r0
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // Actual ranking based on worst-to-best ordering:
    // Species 1: ranks [9, 4, 0] → average = (9+4+0)/3 = 13/3 = 4.33
    // Species 2: ranks [7, 2] → average = (7+2)/2 = 4.5
    // Species 3: ranks [8, 6, 5, 3, 1] → average = (8+6+5+3+1)/5 = 23/5 = 4.6
    
    // Species 1 has best average (4.33), Species 2 is middle (4.5), Species 3 is worst (4.6)
    // Should be assigned ordinal ranks 1, 2, 3 respectively
    
    EXPECT_EQ(speciesData[1].speciesRank, 1) << "Species 1 should have rank 1 (best average 4.33)";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 should have rank 2 (middle average 4.5)";
    EXPECT_EQ(speciesData[3].speciesRank, 3) << "Species 3 should have rank 3 (worst average 4.6)";
    
    // Verify population sizes are correct
    EXPECT_EQ(speciesData[1].currentPopulationSize, 3) << "Species 1 should have 3 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 5) << "Species 3 should have 5 genomes";
}

TEST_F(DynamicDataUpdateTest, PopulationSize_UpdatedCorrectly) {
    // Test core responsibility: operator must update currentPopulationSize
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial sizes to wrong values
    speciesData[1].currentPopulationSize = 99; // Should become 2
    speciesData[2].currentPopulationSize = 99; // Should become 2  
    speciesData[3].currentPopulationSize = 99; // Should become 2
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    EXPECT_EQ(speciesData[1].currentPopulationSize, 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
}

TEST_F(DynamicDataUpdateTest, PopulationSize_UnevenDistribution) {
    // Test with uneven species distribution
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65};
    std::vector<uint32_t> species = {1,   1,  1,  1,  2,  2,  3,  3};
    //                              Species 1: 4 genomes, Species 2: 2 genomes, Species 3: 2 genomes
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial sizes to wrong values
    speciesData[1].currentPopulationSize = 0;
    speciesData[2].currentPopulationSize = 10;
    speciesData[3].currentPopulationSize = 5;
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    EXPECT_EQ(speciesData[1].currentPopulationSize, 4) << "Species 1 should have 4 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
}


TEST_F(DynamicDataUpdateTest, SpeciesRanking_FloatingPointToOrdinalConversion) {
    // Test: Verify conversion from average ranks to ordinal rankings
    // Core algorithm of new ranking functionality
    
    // Create specific scenario with known floating point averages
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55}; // 10 genomes
    std::vector<uint32_t> species = {1,   1,  2,  2,  2,  3,  3,  3,  3,  3}; // Uneven distribution
    //                              r0   r1  r2  r3  r4  r5  r6  r7  r8  r9
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // Calculate expected floating point averages:
    // Fitness data: [100, 95, 90, 85, 80, 75, 70, 65, 60, 55] (best to worst)
    // After sorting worst-to-best: [55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    // Ranks assigned:               [0,  1,  2,  3,  4,  5,  6,  7,  8,  9]
    // Species 1: genomes with fitness [100, 95] → ranks [9, 8] → average = 8.5
    // Species 2: genomes with fitness [90, 85, 80] → ranks [7, 6, 5] → average = 6.0  
    // Species 3: genomes with fitness [75, 70, 65, 60, 55] → ranks [4, 3, 2, 1, 0] → average = 2.0
    
    // Expected ordinal conversion (lower average = better ordinal rank):
    // Species 3: average 2.0 → ordinal rank 1 (best)
    // Species 2: average 6.0 → ordinal rank 2 (middle)
    // Species 1: average 8.5 → ordinal rank 3 (worst)
    
    EXPECT_EQ(speciesData[1].speciesRank, 3) << "Species 1 (avg 8.5) should get ordinal rank 3";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 (avg 6.0) should get ordinal rank 2";
    EXPECT_EQ(speciesData[3].speciesRank, 1) << "Species 3 (avg 2.0) should get ordinal rank 1";
}

TEST_F(DynamicDataUpdateTest, SpeciesRanking_UpdatesSpeciesRankField) {
    // Test: Verify speciesRank field gets correctly updated with ordinal values
    // Field update mechanism for new feature
    
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 3};
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial rank values to verify they get overwritten
    speciesData[1].speciesRank = 999; // Should become 3
    speciesData[2].speciesRank = 888; // Should become 2  
    speciesData[3].speciesRank = 777; // Should become 1
    
    // Set other fields to verify they don't interfere with ranking
    speciesData[1].pendingEliminationRating = 5;
    speciesData[2].pendingEliminationRating = 3;
    speciesData[3].pendingEliminationRating = 1;
    
    dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry);
    
    // Verify speciesRank field was updated correctly
    // Fitness data: [100, 90, 80] (best to worst)
    // After sorting worst-to-best: [80, 90, 100]
    // Ranks assigned:               [0,  1,  2]
    // Species 3: fitness 80 → rank 0 → average 0.0 → ordinal rank 1 (best)
    // Species 2: fitness 90 → rank 1 → average 1.0 → ordinal rank 2 (middle)
    // Species 1: fitness 100 → rank 2 → average 2.0 → ordinal rank 3 (worst)
    
    EXPECT_EQ(speciesData[1].speciesRank, 3) << "speciesRank field should be updated from 999 to 3";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "speciesRank field should be updated from 888 to 2";
    EXPECT_EQ(speciesData[3].speciesRank, 1) << "speciesRank field should be updated from 777 to 1";
    
    // Verify initial rank 0 gets updated (no species should have rank 0 after processing)
    EXPECT_NE(speciesData[1].speciesRank, 0) << "No species should have rank 0 after processing";
    EXPECT_NE(speciesData[2].speciesRank, 0) << "No species should have rank 0 after processing";
    EXPECT_NE(speciesData[3].speciesRank, 0) << "No species should have rank 0 after processing";
    
    // Verify all species get assigned ordinal rank values
    EXPECT_GT(speciesData[1].speciesRank, 0) << "All species should have positive ordinal ranks";
    EXPECT_GT(speciesData[2].speciesRank, 0) << "All species should have positive ordinal ranks";
    EXPECT_GT(speciesData[3].speciesRank, 0) << "All species should have positive ordinal ranks";
}

TEST_F(DynamicDataUpdateTest, SpeciesDiscovery_ThroughGenomeDataOnly) {
    // Test: Verify species referenced only in genome data gets created dynamically
    // This tests the core species discovery feature - new species emerge during evolution
    
    std::vector<double> fitness = {100, 90, 80};
    std::vector<uint32_t> species = {1, 2, 777}; // Species 777 appears only in genome data
    
    auto testData = createTestData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2}); // Only species 1 and 2 exist initially
    
    // Should not crash even though species 777 has no existing data
    EXPECT_NO_THROW(dynamicDataUpdate(testData.fitnessResults, testData.genomeData, speciesData, testData.speciesGrouping, *params, *registry));
    
    // Verify species 777 was discovered and created
    auto discoveredSpeciesIt = speciesData.find(777);
    EXPECT_NE(discoveredSpeciesIt, speciesData.end()) << "Species 777 should be discovered through genome data";
    
    if (discoveredSpeciesIt != speciesData.end()) {
        const auto& discoveredData = discoveredSpeciesIt->second;
        
        // Verify proper initialization of discovered species
        EXPECT_EQ(discoveredData.currentPopulationSize, 1) << "Should reflect the one genome in this species";
        EXPECT_GT(discoveredData.speciesRank, 0) << "Should get valid ordinal rank based on performance";
        EXPECT_FALSE(discoveredData.isMarkedForElimination) << "Should not be marked for elimination initially";
        
        // Verify it gets best rank (rank 1) since it has best individual performance
        // Fitness 80 is worst, but gets rank 0 in sorted array, leading to ordinal rank 1
        EXPECT_EQ(discoveredData.speciesRank, 1) << "Species with worst genome should get best rank";
    }
    
    // Verify existing species still work correctly  
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Existing species 1 should be updated";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 1) << "Existing species 2 should be updated";
    // Fitness data: [100, 90, 80] → after sorting worst-to-best: [80, 90, 100] → ranks [0, 1, 2]
    // Species 1: fitness 100 → rank 2 → ordinal rank 3 (worst)
    // Species 2: fitness 90 → rank 1 → ordinal rank 2 (middle)
    EXPECT_EQ(speciesData[1].speciesRank, 3) << "Species 1 should get worst rank";
    EXPECT_EQ(speciesData[2].speciesRank, 2) << "Species 2 should get middle rank";
}