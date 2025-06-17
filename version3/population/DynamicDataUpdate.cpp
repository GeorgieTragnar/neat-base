#include "DynamicDataUpdate.hpp"
#include "ReproductiveInstruction.hpp"
#include "../logger/Logger.hpp"

namespace Population {

// DynamicDataUpdateParams constructor implementation
DynamicDataUpdateParams::DynamicDataUpdateParams(
    uint32_t maxProtectionLimit,
    uint32_t maxSpeciesProtectionRating,
    double protectedTierPercentage,
    uint32_t worstSpeciesCount,
    uint32_t minActiveSpeciesCount
) : _maxProtectionLimit(maxProtectionLimit),
    _maxSpeciesProtectionRating(maxSpeciesProtectionRating),
    _protectedTierPercentage(protectedTierPercentage),
    _worstSpeciesCount(worstSpeciesCount),
    _minActiveSpeciesCount(minActiveSpeciesCount) {
    
    assert(protectedTierPercentage >= 0.0 && protectedTierPercentage <= 1.0);
    assert(worstSpeciesCount > 0);
    assert(minActiveSpeciesCount > 0);
}

// Internal helper function for Phase 1: Instruction Set Counting and Size Updates
void updateInstructionSetSizes(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationInstructionSets& instructionSets
) {
    // auto logger = LOGGER("population.DynamicDataUpdate");
    // LOG_DEBUG("updateInstructionSetSizes: Processing {} species with instruction sets", instructionSets.size());
    
    // Reset all instruction set sizes to 0 first
    for (auto& [speciesId, data] : speciesData) {
        // LOG_DEBUG("  Resetting species {} instructionSetsSize from {} to 0", speciesId, data.instructionSetsSize);
        data.instructionSetsSize = 0;
    }
    
    // Count instruction sets from GenerationPlanner output
    for (const auto& [speciesId, speciesInstructions] : instructionSets) {
        // LOG_DEBUG("  Processing species {} with {} instruction sets", speciesId, speciesInstructions.size());
        // Early assertion for zero instruction species (catches GenerationPlanner errors)
#ifndef NDEBUG
        // Allow empty instruction sets for species marked for elimination
        // assert(!speciesInstructions.empty() && "GenerationPlanner error: species has zero instructions");
#endif
        
        // Validate instruction sets during counting
#ifndef NDEBUG
        for (const auto& instruction : speciesInstructions) {
            assert(instruction.isValid() && "Malformed instruction detected during counting");
        }
#endif
        
        // Check if species exists in species data
        auto speciesIt = speciesData.find(speciesId);
        if (speciesIt != speciesData.end()) {
            // Update existing species
            // LOG_DEBUG("  Updating existing species {} instructionSetsSize to {}", speciesId, speciesInstructions.size());
            speciesIt->second.instructionSetsSize = static_cast<uint32_t>(speciesInstructions.size());
        } else {
            // Species not found in species data but has instruction sets - this indicates a bug
            // GenerationPlanner should only create instruction sets for existing species
            // auto logger = LOGGER("population.DynamicDataUpdate");
            // LOG_ERROR("GenerationPlanner created instruction sets for non-existent species {}", speciesId);
            assert(false && "GenerationPlanner error: instruction sets created for non-existent species");
        }
    }
    
    // Validate consistency - eliminated species should have zero instruction sets
#ifndef NDEBUG
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination) {
            // Find in instruction sets to verify it has empty or no entry
            auto instructionIt = instructionSets.find(speciesId);
            if (instructionIt != instructionSets.end()) {
                assert(instructionIt->second.empty() && "Eliminated species should have empty instruction sets");
            }
            // If not found in instruction sets, that's also valid (no instructions planned)
        }
    }
#endif
}

// Internal helper function for species ranking based on average performance
void updateSpeciesRanking(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const std::unordered_map<uint32_t, double>& speciesAverageRanks
) {
    // Create vector of (speciesId, averageRank) pairs for sorting
    std::vector<std::pair<uint32_t, double>> speciesRankings;
    for (const auto& [speciesId, avgRank] : speciesAverageRanks) {
        speciesRankings.emplace_back(speciesId, avgRank);
    }
    
    // Sort by average rank (lower rank = better performance)
    std::sort(speciesRankings.begin(), speciesRankings.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    // Assign ordinal rankings (1-based)
    for (size_t i = 0; i < speciesRankings.size(); ++i) {
        uint32_t speciesId = speciesRankings[i].first;
        auto speciesIt = speciesData.find(speciesId);
        if (speciesIt != speciesData.end()) {
            speciesIt->second.speciesRank = static_cast<uint32_t>(i + 1); // 1-based ranking
        }
    }
    
    // Handle species not in rankings (no genomes) - assign worst rank
    uint32_t worstRank = static_cast<uint32_t>(speciesRankings.size() + 1);
    for (auto& [speciesId, data] : speciesData) {
        if (data.speciesRank == 0) { // Not yet assigned
            data.speciesRank = worstRank;
        }
    }
}

} // namespace Population