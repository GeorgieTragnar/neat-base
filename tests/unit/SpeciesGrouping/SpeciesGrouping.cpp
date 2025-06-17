#include <gtest/gtest.h>
#include <map>
#include <unordered_map>

#include "tests/test_common.h"
#include "version3/population/SpeciesGrouping.hpp"

using namespace Population;

class SpeciesGroupingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Test setup if needed
    }
    
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
// 1. BASIC FUNCTIONALITY TESTS
// =============================================================================

TEST_F(SpeciesGroupingTest, BasicGrouping) {
    // Test basic species grouping functionality
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Verify basic structure
    EXPECT_EQ(result.size(), 3) << "Should have exactly 3 species";
    EXPECT_TRUE(result.find(1) != result.end()) << "Species 1 should be present";
    EXPECT_TRUE(result.find(2) != result.end()) << "Species 2 should be present";
    EXPECT_TRUE(result.find(3) != result.end()) << "Species 3 should be present";
    
    // Verify species sizes
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
    
    // Verify indices are in ascending order within species
    EXPECT_EQ(result[1][0], 0) << "Species 1 first genome at index 0";
    EXPECT_EQ(result[1][1], 1) << "Species 1 second genome at index 1";
    EXPECT_EQ(result[2][0], 2) << "Species 2 first genome at index 2";
    EXPECT_EQ(result[2][1], 3) << "Species 2 second genome at index 3";
    EXPECT_EQ(result[3][0], 4) << "Species 3 first genome at index 4";
    EXPECT_EQ(result[3][1], 5) << "Species 3 second genome at index 5";
    
    // Verify population sizes are updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 2) << "Species 1 population size updated";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 population size updated";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 population size updated";
}

TEST_F(SpeciesGroupingTest, SingleSpecies) {
    // Test with all genomes in single species
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {42,  42, 42, 42, 42};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({42});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 1) << "Should have exactly 1 species";
    EXPECT_TRUE(result.find(42) != result.end()) << "Species 42 should be present";
    EXPECT_EQ(result[42].size(), 5) << "Species 42 should have all 5 genomes";
    
    // Verify sequential indices
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ(result[42][i], i) << "Index " << i << " should be at position " << i;
    }
    
    EXPECT_EQ(speciesData[42].currentPopulationSize, 5) << "Population size should be 5";
}

TEST_F(SpeciesGroupingTest, ManySmallSpecies) {
    // Test with many species having few genomes each
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30};
    std::vector<uint32_t> species = {1,   2,  3,  4,  5,  6,  7,  8};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4, 5, 6, 7, 8});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 8) << "Should have exactly 8 species";
    
    // Each species should have exactly 1 genome
    for (uint32_t speciesId = 1; speciesId <= 8; ++speciesId) {
        EXPECT_TRUE(result.find(speciesId) != result.end()) 
            << "Species " << speciesId << " should be present";
        EXPECT_EQ(result[speciesId].size(), 1) 
            << "Species " << speciesId << " should have 1 genome";
        EXPECT_EQ(result[speciesId][0], speciesId - 1) 
            << "Species " << speciesId << " should have genome at index " << (speciesId - 1);
        EXPECT_EQ(speciesData[speciesId].currentPopulationSize, 1) 
            << "Species " << speciesId << " population size should be 1";
    }
}

TEST_F(SpeciesGroupingTest, ImbalancedDistribution) {
    // Test with uneven species sizes
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4};
    //                              Species 1: 8 genomes, Species 2: 3 genomes, Species 3: 2 genomes, Species 4: 3 genomes
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 4) << "Should have exactly 4 species";
    EXPECT_EQ(result[1].size(), 8) << "Species 1 should have 8 genomes";
    EXPECT_EQ(result[2].size(), 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
    EXPECT_EQ(result[4].size(), 3) << "Species 4 should have 3 genomes";
    
    // Verify population sizes
    EXPECT_EQ(speciesData[1].currentPopulationSize, 8);
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3);
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2);
    EXPECT_EQ(speciesData[4].currentPopulationSize, 3);
}

// =============================================================================
// 2. EDGE CASES AND BOUNDARY CONDITIONS
// =============================================================================

TEST_F(SpeciesGroupingTest, SingleGenome) {
    // Test absolute minimum case
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 1) << "Should have exactly 1 species";
    EXPECT_EQ(result[1].size(), 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(result[1][0], 0) << "Single genome should be at index 0";
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1) << "Population size should be 1";
}

TEST_F(SpeciesGroupingTest, ExtremeSpeciesIDs) {
    // Test with boundary species ID values
    std::vector<double> fitness = {100, 90, 80, 70};
    std::vector<uint32_t> species = {0, UINT32_MAX, UINT32_MAX/2, 1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({0, UINT32_MAX, UINT32_MAX/2, 1});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 4) << "Should handle extreme species IDs";
    EXPECT_EQ(result[0].size(), 1) << "Species 0 should work";
    EXPECT_EQ(result[UINT32_MAX].size(), 1) << "Species UINT32_MAX should work";
    EXPECT_EQ(result[UINT32_MAX/2].size(), 1) << "Species UINT32_MAX/2 should work";
    EXPECT_EQ(result[1].size(), 1) << "Species 1 should work";
    
    // Verify population sizes for extreme IDs
    EXPECT_EQ(speciesData[0].currentPopulationSize, 1);
    EXPECT_EQ(speciesData[UINT32_MAX].currentPopulationSize, 1);
    EXPECT_EQ(speciesData[UINT32_MAX/2].currentPopulationSize, 1);
    EXPECT_EQ(speciesData[1].currentPopulationSize, 1);
}

TEST_F(SpeciesGroupingTest, NonSequentialSpeciesIDs) {
    // Test with non-sequential species IDs
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10};
    std::vector<uint32_t> species = {5,   99, 1,  42, 5,  99, 1,  42, 5,  99};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({5, 99, 1, 42});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    EXPECT_EQ(result.size(), 4) << "Should have exactly 4 species";
    EXPECT_EQ(result[5].size(), 3) << "Species 5 should have 3 genomes";
    EXPECT_EQ(result[99].size(), 3) << "Species 99 should have 3 genomes";
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[42].size(), 2) << "Species 42 should have 2 genomes";
}

TEST_F(SpeciesGroupingTest, DuplicateFitnessValues) {
    // Test with duplicate fitness values (multimap behavior)
    std::vector<double> fitness = {100, 100, 100, 90, 90, 80};
    std::vector<uint32_t> species = {1,   2,   3,   1,  2,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Capture the multimap iteration order to verify ordering preservation
    std::vector<uint32_t> multimapSpeciesOrder;
    std::vector<double> multimapFitnessOrder;
    for (const auto& [fitnessResult, genomeData] : genomeData) {
        multimapSpeciesOrder.push_back(genomeData.speciesId);
        multimapFitnessOrder.push_back(fitnessResult.fitness);
    }
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Basic structure validation
    EXPECT_EQ(result.size(), 3) << "Should have exactly 3 species";
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
    
    // Verify ordering preservation: reconstruct the species order from result indices
    std::vector<uint32_t> reconstructedSpeciesOrder(multimapSpeciesOrder.size());
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            reconstructedSpeciesOrder[index] = speciesId;
        }
    }
    
    // The reconstructed order should exactly match the multimap iteration order
    EXPECT_EQ(reconstructedSpeciesOrder, multimapSpeciesOrder) 
        << "Species ordering should be preserved from multimap iteration";
    
    // Verify indices within each species are in ascending order
    for (const auto& [speciesId, indices] : result) {
        for (size_t i = 1; i < indices.size(); ++i) {
            EXPECT_LT(indices[i-1], indices[i]) 
                << "Species " << speciesId << " indices should be in ascending order";
        }
    }
    
    // Verify all indices are used exactly once
    std::set<size_t> allIndices;
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            allIndices.insert(index);
        }
    }
    EXPECT_EQ(allIndices.size(), 6) << "All 6 indices should be present exactly once";
}

// =============================================================================
// 3. DATA INTEGRITY TESTS
// =============================================================================

TEST_F(SpeciesGroupingTest, CompletenessInvariant) {
    // Test systematic verification using set operations that every input genome appears exactly once across all species
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Build complete set of all indices from all species vectors
    std::set<size_t> allIndicesSet;
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            auto [iter, wasInserted] = allIndicesSet.insert(index);
            EXPECT_TRUE(wasInserted) << "Index " << index << " appears multiple times (duplicate detected)";
        }
    }
    
    // Verify set size equals input genome count
    EXPECT_EQ(allIndicesSet.size(), 20) << "Complete index set should contain exactly 20 unique indices";
    EXPECT_EQ(allIndicesSet.size(), fitness.size()) << "Index set size should match input fitness vector size";
    
    // Verify set contents exactly equals {0, 1, 2, ..., n-1}
    std::set<size_t> expectedIndices;
    for (size_t i = 0; i < fitness.size(); ++i) {
        expectedIndices.insert(i);
    }
    EXPECT_EQ(allIndicesSet, expectedIndices) << "Index set should contain exactly {0, 1, 2, ..., 19}";
    
    // Verify no species vector contains duplicate indices within itself
    for (const auto& [speciesId, indices] : result) {
        std::set<size_t> speciesIndicesSet(indices.begin(), indices.end());
        EXPECT_EQ(speciesIndicesSet.size(), indices.size()) 
            << "Species " << speciesId << " contains duplicate indices within its vector";
    }
}

TEST_F(SpeciesGroupingTest, SpeciesIdentityPreservation) {
    // Test bidirectional validation that species assignment is preserved correctly
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45};
    std::vector<uint32_t> species = {1,   2,  3,  4,  1,  2,  3,  4,  1,  2,  3,  4};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Build mapping from global index to original genome data for reference
    std::unordered_map<size_t, uint32_t> indexToOriginalSpecies;
    size_t globalIndex = 0;
    for (const auto& [fitnessResult, genome] : genomeData) {
        indexToOriginalSpecies[globalIndex] = genome.speciesId;
        globalIndex++;
    }
    
    // Forward validation: For each species vector in result, verify all genomes belong to that species
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            EXPECT_TRUE(indexToOriginalSpecies.find(index) != indexToOriginalSpecies.end()) 
                << "Invalid index " << index << " in species " << speciesId;
            EXPECT_EQ(indexToOriginalSpecies[index], speciesId) 
                << "Genome at index " << index << " should belong to species " << speciesId;
        }
    }
    
    // Reverse validation: For each original genome, verify it appears in correct species vector
    for (const auto& [globalIdx, originalSpeciesId] : indexToOriginalSpecies) {
        bool foundInCorrectSpecies = false;
        for (const auto& [speciesId, indices] : result) {
            auto it = std::find(indices.begin(), indices.end(), globalIdx);
            if (it != indices.end()) {
                EXPECT_EQ(speciesId, originalSpeciesId) 
                    << "Genome " << globalIdx << " found in wrong species " << speciesId 
                    << " (should be in " << originalSpeciesId << ")";
                foundInCorrectSpecies = (speciesId == originalSpeciesId);
                break;
            }
        }
        EXPECT_TRUE(foundInCorrectSpecies) 
            << "Genome " << globalIdx << " not found in its correct species " << originalSpeciesId;
    }
}

TEST_F(SpeciesGroupingTest, OrderingPreservation) {
    // Test that fitness ordering is preserved within species
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   1,  2,  2,  1,  2,  3,  3,  3,  1};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Verify indices within each species are in ascending order (fitness rank order)
    for (const auto& [speciesId, indices] : result) {
        for (size_t i = 1; i < indices.size(); ++i) {
            EXPECT_LT(indices[i-1], indices[i]) 
                << "Species " << speciesId << " indices should be in ascending order";
        }
    }
}

// =============================================================================
// 4. TEMPLATE COMPATIBILITY TESTS
// =============================================================================

TEST_F(SpeciesGroupingTest, DifferentFitnessTypes) {
    // Test template compatibility with different fitness result types
    
    // Integer fitness type
    struct IntFitnessResult {
        int fitness;
        bool operator<(const IntFitnessResult& other) const {
            return fitness > other.fitness; // Higher fitness is better
        }
    };
    
    std::multimap<IntFitnessResult, DynamicGenomeData> intGenomeData;
    std::vector<int> intFitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species1 = {1,   1,  2,  2,  3,  3};
    
    for (size_t i = 0; i < intFitness.size(); ++i) {
        IntFitnessResult fitness{intFitness[i]};
        DynamicGenomeData data;
        data.speciesId = species1[i];
        data.protectionCounter = 0;
        data.isUnderRepair = false;
        data.isMarkedForElimination = false;
        
        intGenomeData.insert({fitness, data});
    }
    
    auto speciesData = createTestSpeciesData({1, 2, 3});
    auto result = speciesGrouping(intGenomeData, speciesData);
    
    EXPECT_EQ(result.size(), 3) << "Integer fitness should work";
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
}

// =============================================================================
// 5. POPULATION SIZE UPDATE TESTS
// =============================================================================

TEST_F(SpeciesGroupingTest, PopulationSizeReset) {
    // Test that population sizes are properly reset and updated
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3});
    
    // Set initial population sizes to wrong values
    speciesData[1].currentPopulationSize = 99;
    speciesData[2].currentPopulationSize = 99;
    speciesData[3].currentPopulationSize = 99;
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Verify population sizes are correctly updated
    EXPECT_EQ(speciesData[1].currentPopulationSize, 2) << "Species 1 population should be reset to 2";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 2) << "Species 2 population should be reset to 2";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 population should be reset to 2";
}

TEST_F(SpeciesGroupingTest, PopulationSizeAccuracy) {
    // Test population size accuracy with various distributions
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  4,  4};
    
    auto genomeData = createTestGenomeData(fitness, species);
    auto speciesData = createTestSpeciesData({1, 2, 3, 4});
    
    auto result = speciesGrouping(genomeData, speciesData);
    
    // Verify population sizes match actual genome counts
    EXPECT_EQ(speciesData[1].currentPopulationSize, 5) << "Species 1 should have 5 genomes";
    EXPECT_EQ(speciesData[2].currentPopulationSize, 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
    EXPECT_EQ(speciesData[4].currentPopulationSize, 5) << "Species 4 should have 5 genomes";
    
    // Cross-verify with result indices
    EXPECT_EQ(result[1].size(), speciesData[1].currentPopulationSize);
    EXPECT_EQ(result[2].size(), speciesData[2].currentPopulationSize);
    EXPECT_EQ(result[3].size(), speciesData[3].currentPopulationSize);
    EXPECT_EQ(result[4].size(), speciesData[4].currentPopulationSize);
}