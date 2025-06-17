#pragma once

#include <unordered_map>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdint>

#include "PopulationData.hpp"
#include "GenerationPlannerParams.hpp"
#include "ReproductiveInstruction.hpp"

namespace Population {

// Population health assessment results from Phase 1
struct PopulationHealthAssessment {
    uint32_t activeSpeciesCount;        // Number of non-eliminated species
    uint32_t totalCurrentInstructions; // Sum of all instructionSetsSize values
    uint32_t equilibriumTarget;         // Target instruction sets per species
    float averageSpeciesSize;           // Mean species instruction set size
    float sizeVariance;                 // Variance in species sizes
    bool isSystemHealthy;               // Overall system health indicator
};

// Elite allocation results from Phase 2
struct EliteAllocationResult {
    std::unordered_map<uint32_t, uint32_t> eliteCountPerSpecies; // speciesId -> elite count
    uint32_t totalEliteCount = 0;                                // Total elites across all species
};

// Crossover allocation results from Phase 3
struct CrossoverAllocationResult {
    std::unordered_map<uint32_t, uint32_t> crossoverCountPerSpecies; // speciesId -> crossover count
    uint32_t totalCrossoverCount = 0;                                 // Total crossovers across all species
};

// Equilibrium-driven mutation allocation results from Phase 4
struct EquilibriumMutationAllocationResult {
    std::unordered_map<uint32_t, uint32_t> guaranteedUnprotectedMutationsPerSpecies; // speciesId -> guaranteed unprotected count
    uint32_t totalGuaranteedUnprotectedMutations = 0;                                 // Total guaranteed unprotected across all species
};

// Protection assignment results from Phase 5
struct ProtectionAssignmentResult {
    std::unordered_map<uint32_t, uint32_t> unprotectedRemainderPerSpecies; // speciesId -> unprotected remainder count
    std::unordered_map<uint32_t, uint32_t> protectedRemainderPerSpecies;   // speciesId -> protected remainder count
    uint32_t totalUnprotectedRemainder = 0;                                 // Total unprotected remainder across all species
    uint32_t totalProtectedRemainder = 0;                                   // Total protected remainder across all species
};

// Function declarations

// Phase 1: Population Health Assessment
PopulationHealthAssessment assessPopulationHealth(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
);

// Phase 2: Elite Preservation Allocation with rank-based scaling
EliteAllocationResult allocateElitePreservation(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
);

// Phase 3: Crossover Allocation with constraints and scaling
CrossoverAllocationResult allocateCrossover(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EliteAllocationResult& eliteAllocation,
    const GenerationPlannerParams& params
);

// Phase 4: Equilibrium-Driven Mutation Allocation
EquilibriumMutationAllocationResult allocateEquilibriumMutations(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const PopulationHealthAssessment& healthAssessment,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const GenerationPlannerParams& params
);

// Phase 5: Protection Assignment for Remaining Genomes
ProtectionAssignmentResult assignProtectionForRemainingGenomes(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const EquilibriumMutationAllocationResult& equilibriumMutationAllocation,
    const GenerationPlannerParams& params
);

// Phase 6: Instruction Set Generation with exact count matching
GenerationInstructionSets generateInstructionSets(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const EquilibriumMutationAllocationResult& equilibriumMutationAllocation,
    const ProtectionAssignmentResult& protectionAssignment,
    const GenerationPlannerParams& params
);

// Main GenerationPlanner operator function
GenerationInstructionSets generationPlanner(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
);

} // namespace Population