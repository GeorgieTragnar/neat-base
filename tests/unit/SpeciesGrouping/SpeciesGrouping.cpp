#include <gtest/gtest.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "tests/test_common.h"
#include "version3/population/SpeciesGrouping.hpp"
#include "version3/population/GlobalIndexRegistry.hpp"

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
    
    // Helper to create test data matching the current function signature
    struct TestData {
        std::multimap<TestFitnessResult, size_t> fitnessResults;
        std::vector<DynamicGenomeData> genomeData;
        std::unordered_map<uint32_t, DynamicSpeciesData> speciesData;
        GlobalIndexRegistry registry;
        
        TestData(uint32_t maxGenomes) : registry(maxGenomes) {}
    };
    
    // Helper to create complete test data set
    TestData createTestData(
        const std::vector<double>& fitnessValues,
        const std::vector<uint32_t>& speciesIds,
        const std::vector<GenomeState>& genomeStates = {},
        const std::vector<bool>& underRepair = {},
        const std::vector<uint32_t>& eliminatedSpecies = {}) {
        
        EXPECT_EQ(fitnessValues.size(), speciesIds.size());
        
        TestData data(static_cast<uint32_t>(fitnessValues.size()));
        
        // Create genome data vector
        data.genomeData.resize(fitnessValues.size());
        for (size_t i = 0; i < fitnessValues.size(); ++i) {
            data.genomeData[i].speciesId = speciesIds[i];
            data.genomeData[i].isUnderRepair = underRepair.empty() ? false : underRepair[i];
            data.genomeData[i].isMarkedForElimination = false;
            data.genomeData[i].pendingEliminationCounter = 0;
            data.genomeData[i].repairAttempts = 0;
            data.genomeData[i].genomeIndex = static_cast<uint32_t>(i);
            data.genomeData[i].parentAIndex = UINT32_MAX;
            data.genomeData[i].parentBIndex = UINT32_MAX;
        }
        
        // Create fitness results multimap
        for (size_t i = 0; i < fitnessValues.size(); ++i) {
            TestFitnessResult fitness{fitnessValues[i]};
            data.fitnessResults.insert({fitness, i});
        }
        
        // Set up registry states
        for (size_t i = 0; i < fitnessValues.size(); ++i) {
            GenomeState state = genomeStates.empty() ? GenomeState::Active : genomeStates[i];
            if (state == GenomeState::Elite) {
                data.registry.markAsElite(static_cast<uint32_t>(i));
            } else if (state == GenomeState::HotElimination) {
                data.registry.markForElimination(static_cast<uint32_t>(i));
            } else if (state == GenomeState::ColdElimination) {
                // Transition from Active -> HotElimination -> ColdElimination
                data.registry.markForElimination(static_cast<uint32_t>(i));
                data.registry.transitionToCold(static_cast<uint32_t>(i));
            } else if (state == GenomeState::ReadyForReplacement) {
                // Transition from Active -> HotElimination -> ColdElimination -> ReadyForReplacement
                data.registry.markForElimination(static_cast<uint32_t>(i));
                data.registry.transitionToCold(static_cast<uint32_t>(i));
                data.registry.markReadyForReplacement(static_cast<uint32_t>(i));
            }
        }
        
        // Create species data
        std::set<uint32_t> allSpeciesIds(speciesIds.begin(), speciesIds.end());
        for (uint32_t id : allSpeciesIds) {
            DynamicSpeciesData speciesInfo;
            speciesInfo.pendingEliminationRating = 0;
            speciesInfo.currentPopulationSize = 0;
            speciesInfo.speciesRank = 0;
            speciesInfo.isMarkedForElimination = 
                std::find(eliminatedSpecies.begin(), eliminatedSpecies.end(), id) != eliminatedSpecies.end();
            
            data.speciesData[id] = speciesInfo;
        }
        
        return data;
    }
};

TEST_F(SpeciesGroupingTest, BasicGrouping) {
    // Test basic species grouping functionality with all valid genomes
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Verify basic structure
    EXPECT_EQ(result.size(), 3) << "Should have exactly 3 species";
    EXPECT_TRUE(result.find(1) != result.end()) << "Species 1 should be present";
    EXPECT_TRUE(result.find(2) != result.end()) << "Species 2 should be present";
    EXPECT_TRUE(result.find(3) != result.end()) << "Species 3 should be present";
    
    // Verify species sizes
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
    
    // Verify fitness ordering is preserved (higher fitness comes first)
    EXPECT_LT(result[1][0], result[1][1]) << "Species 1 should be ordered by fitness rank";
    EXPECT_LT(result[2][0], result[2][1]) << "Species 2 should be ordered by fitness rank";
    EXPECT_LT(result[3][0], result[3][1]) << "Species 3 should be ordered by fitness rank";
    
    // Verify population sizes are updated
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 2) << "Species 1 population size updated";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 2) << "Species 2 population size updated";
    EXPECT_EQ(data.speciesData[3].currentPopulationSize, 2) << "Species 3 population size updated";
}

TEST_F(SpeciesGroupingTest, SingleSpecies) {
    // Test with all genomes in single species
    std::vector<double> fitness = {100, 90, 80, 70, 60};
    std::vector<uint32_t> species = {42,  42, 42, 42, 42};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    EXPECT_EQ(result.size(), 1) << "Should have exactly 1 species";
    EXPECT_TRUE(result.find(42) != result.end()) << "Species 42 should be present";
    EXPECT_EQ(result[42].size(), 5) << "Species 42 should have all 5 genomes";
    
    // Verify fitness ordering is preserved
    for (size_t i = 1; i < result[42].size(); ++i) {
        EXPECT_LT(result[42][i-1], result[42][i]) << "Genomes should be ordered by fitness rank";
    }
    
    EXPECT_EQ(data.speciesData[42].currentPopulationSize, 5) << "Population size should be 5";
}

TEST_F(SpeciesGroupingTest, SingleGenome) {
    // Test absolute minimum case
    std::vector<double> fitness = {100};
    std::vector<uint32_t> species = {1};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    EXPECT_EQ(result.size(), 1) << "Should have exactly 1 species";
    EXPECT_EQ(result[1].size(), 1) << "Species 1 should have 1 genome";
    EXPECT_EQ(result[1][0], 0) << "Single genome should be at index 0";
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 1) << "Population size should be 1";
}

TEST_F(SpeciesGroupingTest, NonSequentialSpeciesIDs) {
    // Test with non-sequential species IDs
    std::vector<double> fitness = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10};
    std::vector<uint32_t> species = {5,   99, 1,  42, 5,  99, 1,  42, 5,  99};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    EXPECT_EQ(result.size(), 4) << "Should have exactly 4 species";
    EXPECT_EQ(result[5].size(), 3) << "Species 5 should have 3 genomes";
    EXPECT_EQ(result[99].size(), 3) << "Species 99 should have 3 genomes";
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[42].size(), 2) << "Species 42 should have 2 genomes";
}

TEST_F(SpeciesGroupingTest, CompletenessInvariant) {
    // Test that every valid genome appears exactly once across all species
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10, 5};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Build complete set of all indices from all species vectors
    std::set<size_t> allIndicesSet;
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            auto [iter, wasInserted] = allIndicesSet.insert(index);
            EXPECT_TRUE(wasInserted) << "Index " << index << " appears multiple times (duplicate detected)";
        }
    }
    
    // All genomes should be valid (Active state, not under repair, not from eliminated species)
    // so completeness means all genomes appear exactly once
    EXPECT_EQ(allIndicesSet.size(), fitness.size()) << "All valid genomes should appear exactly once";
    
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
    // Test that species assignment is preserved correctly
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45};
    std::vector<uint32_t> species = {1,   2,  3,  4,  1,  2,  3,  4,  1,  2,  3,  4};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Forward validation: For each species vector in result, verify all genomes belong to that species
    for (const auto& [speciesId, indices] : result) {
        for (size_t index : indices) {
            EXPECT_LT(index, data.genomeData.size()) << "Invalid index " << index;
            EXPECT_EQ(data.genomeData[index].speciesId, speciesId) 
                << "Genome at index " << index << " should belong to species " << speciesId;
        }
    }
    
    // Reverse validation: For each valid genome, verify it appears in correct species vector
    for (size_t i = 0; i < data.genomeData.size(); ++i) {
        uint32_t expectedSpeciesId = data.genomeData[i].speciesId;
        
        bool foundInCorrectSpecies = false;
        for (const auto& [speciesId, indices] : result) {
            auto it = std::find(indices.begin(), indices.end(), i);
            if (it != indices.end()) {
                EXPECT_EQ(speciesId, expectedSpeciesId) 
                    << "Genome " << i << " found in wrong species " << speciesId 
                    << " (should be in " << expectedSpeciesId << ")";
                foundInCorrectSpecies = (speciesId == expectedSpeciesId);
                break;
            }
        }
        EXPECT_TRUE(foundInCorrectSpecies) 
            << "Genome " << i << " not found in its correct species " << expectedSpeciesId;
    }
}

TEST_F(SpeciesGroupingTest, FitnessOrderingPreservation) {
    // Test that fitness ordering is preserved within species
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55};
    std::vector<uint32_t> species = {1,   1,  2,  2,  1,  2,  3,  3,  3,  1};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Verify indices within each species are in fitness rank order (ascending index = higher fitness)
    for (const auto& [speciesId, indices] : result) {
        for (size_t i = 1; i < indices.size(); ++i) {
            EXPECT_LT(indices[i-1], indices[i]) 
                << "Species " << speciesId << " indices should preserve fitness ordering";
        }
    }
}

TEST_F(SpeciesGroupingTest, PopulationSizeAccuracy) {
    // Test population size accuracy with various distributions
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30};
    std::vector<uint32_t> species = {1,   1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  4,  4};
    
    auto data = createTestData(fitness, species);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Verify population sizes match actual genome counts
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 5) << "Species 1 should have 5 genomes";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 3) << "Species 2 should have 3 genomes";
    EXPECT_EQ(data.speciesData[3].currentPopulationSize, 2) << "Species 3 should have 2 genomes";
    EXPECT_EQ(data.speciesData[4].currentPopulationSize, 5) << "Species 4 should have 5 genomes";
    
    // Cross-verify with result indices
    EXPECT_EQ(result[1].size(), data.speciesData[1].currentPopulationSize);
    EXPECT_EQ(result[2].size(), data.speciesData[2].currentPopulationSize);
    EXPECT_EQ(result[3].size(), data.speciesData[3].currentPopulationSize);
    EXPECT_EQ(result[4].size(), data.speciesData[4].currentPopulationSize);
}

TEST_F(SpeciesGroupingTest, GenomeStateFiltering) {
    // Test that only Active and Elite genomes are included in result
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  1,  2,  2,  2};
    std::vector<GenomeState> states = {
        GenomeState::Active,        // Index 0 - should be included
        GenomeState::Elite,         // Index 1 - should be included
        GenomeState::HotElimination,// Index 2 - should NOT be included
        GenomeState::Active,        // Index 3 - should be included
        GenomeState::ColdElimination,// Index 4 - should NOT be included
        GenomeState::ReadyForReplacement// Index 5 - should NOT be included
    };
    
    auto data = createTestData(fitness, species, states);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Only Active and Elite genomes should be in result
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 valid genomes (Active+Elite)";
    EXPECT_EQ(result[2].size(), 1) << "Species 2 should have 1 valid genome (Active only)";
    
    // Verify specific indices
    EXPECT_EQ(result[1][0], 0) << "Species 1 should contain Active genome at index 0";
    EXPECT_EQ(result[1][1], 1) << "Species 1 should contain Elite genome at index 1";
    EXPECT_EQ(result[2][0], 3) << "Species 2 should contain Active genome at index 3";
    
    // Population counts should include ALL genomes, not just valid ones
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 3) << "Species 1 total population includes invalid genomes";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 3) << "Species 2 total population includes invalid genomes";
}

TEST_F(SpeciesGroupingTest, UnderRepairFiltering) {
    // Test that genomes under repair are excluded even if they have valid states
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  1,  2,  2,  2};
    std::vector<bool> underRepair = {false, true, false, true, false, false};
    
    auto data = createTestData(fitness, species, {}, underRepair);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Only non-repair genomes should be in result
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 non-repair genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 non-repair genomes";
    
    // Verify specific indices (should skip repair genomes)
    EXPECT_EQ(result[1][0], 0) << "Species 1 should contain non-repair genome at index 0";
    EXPECT_EQ(result[1][1], 2) << "Species 1 should contain non-repair genome at index 2";
    EXPECT_EQ(result[2][0], 4) << "Species 2 should contain non-repair genome at index 4";
    EXPECT_EQ(result[2][1], 5) << "Species 2 should contain non-repair genome at index 5";
    
    // Population counts should include ALL genomes
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 3) << "Species 1 total population includes repair genomes";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 3) << "Species 2 total population includes repair genomes";
}

TEST_F(SpeciesGroupingTest, EliminatedSpeciesFiltering) {
    // Test that genomes from eliminated species are excluded from result
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  2,  2,  3,  3};
    std::vector<uint32_t> eliminatedSpecies = {2}; // Species 2 is marked for elimination
    
    auto data = createTestData(fitness, species, {}, {}, eliminatedSpecies);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Only non-eliminated species should be in result
    EXPECT_EQ(result.size(), 2) << "Should have 2 non-eliminated species";
    EXPECT_TRUE(result.find(1) != result.end()) << "Species 1 should be present";
    EXPECT_TRUE(result.find(2) == result.end()) << "Species 2 should NOT be present (eliminated)";
    EXPECT_TRUE(result.find(3) != result.end()) << "Species 3 should be present";
    
    // Verify species sizes
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[3].size(), 2) << "Species 3 should have 2 genomes";
    
    // Population counts should still include ALL genomes, even from eliminated species
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 2) << "Species 1 population count";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 2) << "Species 2 population count (eliminated but still counted)";
    EXPECT_EQ(data.speciesData[3].currentPopulationSize, 2) << "Species 3 population count";
}

TEST_F(SpeciesGroupingTest, NewSpeciesDiscovery) {
    // Test that genomes with species IDs not in speciesData are still grouped
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75};
    std::vector<uint32_t> species = {1,   1,  2,  2,  99, 99}; // Species 99 is new
    
    auto data = createTestData(fitness, species);
    // Remove species 99 from speciesData to simulate new species discovery
    data.speciesData.erase(99);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // All species should be in result, including the new one
    EXPECT_EQ(result.size(), 3) << "Should have 3 species including new discovery";
    EXPECT_TRUE(result.find(1) != result.end()) << "Species 1 should be present";
    EXPECT_TRUE(result.find(2) != result.end()) << "Species 2 should be present";
    EXPECT_TRUE(result.find(99) != result.end()) << "Species 99 should be present (new discovery)";
    
    // Verify species sizes
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 genomes";
    EXPECT_EQ(result[2].size(), 2) << "Species 2 should have 2 genomes";
    EXPECT_EQ(result[99].size(), 2) << "Species 99 should have 2 genomes";
    
    // Population counts should only be updated for existing species
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 2) << "Species 1 population count";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 2) << "Species 2 population count";
    EXPECT_TRUE(data.speciesData.find(99) == data.speciesData.end()) << "Species 99 should not be in speciesData";
}

TEST_F(SpeciesGroupingTest, CombinedFilteringScenario) {
    // Test comprehensive filtering: states, repair, and elimination combined
    std::vector<double> fitness = {100, 95, 90, 85, 80, 75, 70, 65};
    std::vector<uint32_t> species = {1,   1,  1,  1,  2,  2,  2,  2};
    std::vector<GenomeState> states = {
        GenomeState::Active,        // Index 0 - valid
        GenomeState::Elite,         // Index 1 - valid
        GenomeState::HotElimination,// Index 2 - invalid state
        GenomeState::Active,        // Index 3 - valid
        GenomeState::Active,        // Index 4 - valid
        GenomeState::Elite,         // Index 5 - valid
        GenomeState::Active,        // Index 6 - valid
        GenomeState::Active         // Index 7 - valid
    };
    std::vector<bool> underRepair = {false, false, false, true, false, false, true, false};
    std::vector<uint32_t> eliminatedSpecies = {2}; // Species 2 is eliminated
    
    auto data = createTestData(fitness, species, states, underRepair, eliminatedSpecies);
    
    auto result = speciesGrouping(data.fitnessResults, data.genomeData, data.speciesData, data.registry);
    
    // Only species 1 should be in result (species 2 is eliminated)
    EXPECT_EQ(result.size(), 1) << "Should have 1 non-eliminated species";
    EXPECT_TRUE(result.find(1) != result.end()) << "Species 1 should be present";
    EXPECT_TRUE(result.find(2) == result.end()) << "Species 2 should NOT be present (eliminated)";
    
    // Species 1 should have only valid genomes: indices 0, 1 (index 2 invalid state, index 3 under repair)
    EXPECT_EQ(result[1].size(), 2) << "Species 1 should have 2 valid genomes";
    EXPECT_EQ(result[1][0], 0) << "Species 1 should contain genome at index 0";
    EXPECT_EQ(result[1][1], 1) << "Species 1 should contain genome at index 1";
    
    // Population counts should include ALL genomes regardless of validity
    EXPECT_EQ(data.speciesData[1].currentPopulationSize, 4) << "Species 1 total population";
    EXPECT_EQ(data.speciesData[2].currentPopulationSize, 4) << "Species 2 total population";
}