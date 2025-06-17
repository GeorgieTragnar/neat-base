#include "GenerationPlanner.hpp"
#include "../logger/Logger.hpp"

namespace Population {

// Phase 1: Population Health Assessment
PopulationHealthAssessment assessPopulationHealth(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
) {
    auto logger = LOGGER("population.GenerationPlanner");
    LOG_DEBUG("assessPopulationHealth: Processing {} species", speciesData.size());
    
    PopulationHealthAssessment assessment = {};
    
    // Count active species and total instruction sets
    std::vector<uint32_t> activeSizes;
    activeSizes.reserve(speciesData.size());
    
    for (const auto& [speciesId, data] : speciesData) {
        LOG_DEBUG("  Species {}: instructionSetsSize={}, isMarkedForElimination={}", 
                 speciesId, data.instructionSetsSize, data.isMarkedForElimination);
        if (!data.isMarkedForElimination) {
            assessment.activeSpeciesCount++;
            assessment.totalCurrentInstructions += data.instructionSetsSize;
            activeSizes.push_back(data.instructionSetsSize);
        }
    }
    
    LOG_DEBUG("  Active species count: {}", assessment.activeSpeciesCount);
    LOG_DEBUG("  Total current instructions: {}", assessment.totalCurrentInstructions);
    
    // Handle edge cases
    if (assessment.activeSpeciesCount == 0) {
        // Emergency fallback - no active species
        assessment.equilibriumTarget = params.getMinSpeciesSize();
        assessment.averageSpeciesSize = 0.0f;
        assessment.sizeVariance = 0.0f;
        assessment.isSystemHealthy = false;
        return assessment;
    }
    
    if (assessment.activeSpeciesCount == 1) {
        // Single species - still apply equilibrium force to reach target population
        uint32_t rawEquilibriumTarget = params.getTargetTotalPopulation(); // All population goes to single species
        assessment.equilibriumTarget = std::max(rawEquilibriumTarget, params.getMinSpeciesSize());
        assessment.averageSpeciesSize = static_cast<float>(activeSizes[0]);
        assessment.sizeVariance = 0.0f;
        assessment.isSystemHealthy = true;
        
        LOG_DEBUG("  Single species mode: Target population: {}, Final equilibrium target: {}", 
                 params.getTargetTotalPopulation(), assessment.equilibriumTarget);
        return assessment;
    }
    
    // Calculate equilibrium target
    uint32_t rawEquilibriumTarget = params.getTargetTotalPopulation() / assessment.activeSpeciesCount;
    assessment.equilibriumTarget = std::max(rawEquilibriumTarget, params.getMinSpeciesSize());
    
    LOG_DEBUG("  Target population: {}, Raw equilibrium target: {}, Min species size: {}", 
             params.getTargetTotalPopulation(), rawEquilibriumTarget, params.getMinSpeciesSize());
    LOG_DEBUG("  Final equilibrium target: {}", assessment.equilibriumTarget);
    
    // Calculate system health metrics
    if (!activeSizes.empty()) {
        // Average species size
        uint32_t totalSize = std::accumulate(activeSizes.begin(), activeSizes.end(), 0u);
        assessment.averageSpeciesSize = static_cast<float>(totalSize) / activeSizes.size();
        
        // Size distribution variance
        float sumSquaredDiffs = 0.0f;
        for (uint32_t size : activeSizes) {
            float diff = static_cast<float>(size) - assessment.averageSpeciesSize;
            sumSquaredDiffs += diff * diff;
        }
        assessment.sizeVariance = sumSquaredDiffs / activeSizes.size();
        
        // System health assessment
        // Healthy if:
        // 1. Average size is reasonably close to equilibrium target (within 50%)
        // 2. Variance is not extremely high (coefficient of variation < 1.0)
        float avgToTargetRatio = assessment.averageSpeciesSize / assessment.equilibriumTarget;
        float coefficientOfVariation = (assessment.averageSpeciesSize > 0.0f) ? 
            std::sqrt(assessment.sizeVariance) / assessment.averageSpeciesSize : 0.0f;
        
        assessment.isSystemHealthy = (avgToTargetRatio >= 0.5f && avgToTargetRatio <= 2.0f) &&
                                   (coefficientOfVariation < 1.0f);
        
        LOG_DEBUG("  Average species size: {:.2f}, Size variance: {:.2f}", 
                 assessment.averageSpeciesSize, assessment.sizeVariance);
        LOG_DEBUG("  Avg/Target ratio: {:.2f}, Coeff of variation: {:.2f}, System healthy: {}", 
                 avgToTargetRatio, coefficientOfVariation, assessment.isSystemHealthy);
    } else {
        assessment.averageSpeciesSize = 0.0f;
        assessment.sizeVariance = 0.0f;
        assessment.isSystemHealthy = false;
    }
    
    #ifdef DEBUG
    // Debug validation
    assert(assessment.activeSpeciesCount <= speciesData.size());
    assert(assessment.equilibriumTarget >= params.getMinSpeciesSize());
    assert(assessment.averageSpeciesSize >= 0.0f);
    assert(assessment.sizeVariance >= 0.0f);
    
    // Validate minimum active species count constraint
    // Note: This may be violated temporarily during population initialization,
    // but should be maintained after the first few generations
    if (speciesData.size() >= params.getMinActiveSpeciesCount()) {
        // Only assert if we have enough total species to meet the minimum
        assert(assessment.activeSpeciesCount >= params.getMinActiveSpeciesCount() && 
               "Active species count is below minimum - check DynamicDataUpdate parameters");
    }
    #endif
    
    return assessment;
}

// Phase 2: Elite Preservation Allocation with rank-based scaling
EliteAllocationResult allocateElitePreservation(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
) {
    EliteAllocationResult result;
    
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination || data.instructionSetsSize == 0) {
            // Eliminated species or empty species get no elites
            result.eliteCountPerSpecies[speciesId] = 0;
            continue;
        }
        
        // Calculate base elite count with rank-based scaling
        float scalingFactor = params.getEliteScalingFactor(data.speciesRank);
        uint32_t scaledEliteCount = static_cast<uint32_t>(
            std::ceil(params.getBaseEliteCount() * scalingFactor)
        );
        
        // Apply maximum elite percentage constraint
        uint32_t maxElitesByPercentage = static_cast<uint32_t>(
            std::floor(data.instructionSetsSize * params.getMaxElitePercentage())
        );
        
        // Final elite count is minimum of scaled count, percentage cap, and species size
        uint32_t eliteCount = std::min({
            scaledEliteCount,
            maxElitesByPercentage,
            data.instructionSetsSize
        });
        
        // Guarantee at least 1 elite per non-empty species (if instructionSetsSize >= 1)
        if (data.instructionSetsSize >= 1 && eliteCount == 0) {
            eliteCount = 1;
        }
        
        result.eliteCountPerSpecies[speciesId] = eliteCount;
        result.totalEliteCount += eliteCount;
        
        #ifdef DEBUG
        // Validation
        assert(eliteCount <= data.instructionSetsSize);
        assert(eliteCount == 0 || data.instructionSetsSize > 0);
        if (data.instructionSetsSize > 0) {
            float actualPercentage = static_cast<float>(eliteCount) / data.instructionSetsSize;
            assert(actualPercentage <= params.getMaxElitePercentage() + 0.001f); // Allow small floating point error
        }
        #endif
    }
    
    return result;
}

// Phase 3: Crossover Allocation with constraints and scaling
CrossoverAllocationResult allocateCrossover(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EliteAllocationResult& eliteAllocation,
    const GenerationPlannerParams& params
) {
    CrossoverAllocationResult result;
    
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination || data.instructionSetsSize == 0) {
            // Eliminated species or empty species get no crossover
            result.crossoverCountPerSpecies[speciesId] = 0;
            continue;
        }
        
        // Check minimum species size requirement for crossover
        if (!params.canEnableCrossover(data.instructionSetsSize)) {
            // Species too small for crossover (requires at least 2 genomes)
            result.crossoverCountPerSpecies[speciesId] = 0;
            continue;
        }
        
        // Calculate base crossover count with rank-based scaling
        float scalingFactor = params.getCrossoverScalingFactor(data.speciesRank);
        uint32_t scaledCrossoverCount = static_cast<uint32_t>(
            std::ceil(params.getBaseCrossoverSlots() * scalingFactor)
        );
        
        // Get elite count for this species to calculate remaining slots
        auto eliteIt = eliteAllocation.eliteCountPerSpecies.find(speciesId);
        uint32_t eliteCount = (eliteIt != eliteAllocation.eliteCountPerSpecies.end()) ? 
                             eliteIt->second : 0;
        
        // Calculate available slots (cannot consume elite-designated slots)
        uint32_t availableSlots = (data.instructionSetsSize > eliteCount) ? 
                                 (data.instructionSetsSize - eliteCount) : 0;
        
        // Final crossover count is minimum of scaled count and available slots
        uint32_t crossoverCount = std::min(scaledCrossoverCount, availableSlots);
        
        result.crossoverCountPerSpecies[speciesId] = crossoverCount;
        result.totalCrossoverCount += crossoverCount;
        
        #ifdef DEBUG
        // Validation
        assert(crossoverCount <= availableSlots);
        assert(eliteCount + crossoverCount <= data.instructionSetsSize);
        if (crossoverCount > 0) {
            assert(data.instructionSetsSize >= params.getMinSpeciesSizeForCrossover());
        }
        #endif
    }
    
    return result;
}

// Phase 4: Equilibrium-Driven Mutation Allocation
EquilibriumMutationAllocationResult allocateEquilibriumMutations(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const PopulationHealthAssessment& healthAssessment,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const GenerationPlannerParams& params
) {
    auto logger = LOGGER("population.GenerationPlanner");
    LOG_DEBUG("allocateEquilibriumMutations: Equilibrium target = {}", healthAssessment.equilibriumTarget);
    
    EquilibriumMutationAllocationResult result;
    
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination || data.instructionSetsSize == 0) {
            // Eliminated species or empty species get no guaranteed unprotected mutations
            result.guaranteedUnprotectedMutationsPerSpecies[speciesId] = 0;
            continue;
        }
        
        // Get current allocation for this species
        auto eliteIt = eliteAllocation.eliteCountPerSpecies.find(speciesId);
        auto crossoverIt = crossoverAllocation.crossoverCountPerSpecies.find(speciesId);
        
        uint32_t eliteCount = (eliteIt != eliteAllocation.eliteCountPerSpecies.end()) ? 
                             eliteIt->second : 0;
        uint32_t crossoverCount = (crossoverIt != crossoverAllocation.crossoverCountPerSpecies.end()) ? 
                                 crossoverIt->second : 0;
        
        uint32_t currentAllocation = eliteCount + crossoverCount;
        
        // Calculate equilibrium gap
        // Positive gap = species below equilibrium (needs growth)
        // Negative/zero gap = species at/above equilibrium (no guaranteed unprotected mutations)
        int32_t equilibriumGap = static_cast<int32_t>(healthAssessment.equilibriumTarget) - 
                                static_cast<int32_t>(currentAllocation);
        
        LOG_DEBUG("  Species {}: currentAllocation={}, equilibriumTarget={}, equilibriumGap={}", 
                 speciesId, currentAllocation, healthAssessment.equilibriumTarget, equilibriumGap);
        
        uint32_t guaranteedUnprotectedMutations = 0;
        
        if (equilibriumGap > 0) {
            // Species below equilibrium - apply equilibrium bias
            float biasedGap = equilibriumGap * params.getEquilibriumBias();
            uint32_t maxPossibleGap = static_cast<uint32_t>(biasedGap);
            
            // For growth scenarios, we allow expansion beyond current instructionSetsSize
            // The limit is the equilibrium target itself
            uint32_t maxAllowedTotal = healthAssessment.equilibriumTarget;
            uint32_t remainingSlots = (maxAllowedTotal > currentAllocation) ? 
                                     (maxAllowedTotal - currentAllocation) : 0;
            
            // Guaranteed unprotected mutations is minimum of biased gap and remaining slots
            guaranteedUnprotectedMutations = std::min(maxPossibleGap, remainingSlots);
            
            LOG_DEBUG("    Equilibrium bias: {:.2f}, biasedGap: {:.2f}, maxPossibleGap: {}, remainingSlots: {}, guaranteedUnprotected: {}", 
                     params.getEquilibriumBias(), biasedGap, maxPossibleGap, remainingSlots, guaranteedUnprotectedMutations);
        }
        // else: equilibriumGap <= 0, so guaranteedUnprotectedMutations remains 0
        
        result.guaranteedUnprotectedMutationsPerSpecies[speciesId] = guaranteedUnprotectedMutations;
        result.totalGuaranteedUnprotectedMutations += guaranteedUnprotectedMutations;
        
        #ifdef DEBUG
        // Validation
        assert(currentAllocation + guaranteedUnprotectedMutations <= data.instructionSetsSize);
        if (guaranteedUnprotectedMutations > 0) {
            assert(equilibriumGap > 0); // Should only get guaranteed mutations when below equilibrium
        }
        #endif
    }
    
    return result;
}

// Phase 5: Protection Assignment for Remaining Genomes
ProtectionAssignmentResult assignProtectionForRemainingGenomes(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const PopulationHealthAssessment& healthAssessment,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const EquilibriumMutationAllocationResult& equilibriumMutationAllocation,
    const GenerationPlannerParams& params
) {
    ProtectionAssignmentResult result;
    
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination || data.instructionSetsSize == 0) {
            // Eliminated species or empty species get no remainder assignments
            result.unprotectedRemainderPerSpecies[speciesId] = 0;
            result.protectedRemainderPerSpecies[speciesId] = 0;
            continue;
        }
        
        // Get current allocations for this species
        auto eliteIt = eliteAllocation.eliteCountPerSpecies.find(speciesId);
        auto crossoverIt = crossoverAllocation.crossoverCountPerSpecies.find(speciesId);
        auto guaranteedIt = equilibriumMutationAllocation.guaranteedUnprotectedMutationsPerSpecies.find(speciesId);
        
        uint32_t eliteCount = (eliteIt != eliteAllocation.eliteCountPerSpecies.end()) ? 
                             eliteIt->second : 0;
        uint32_t crossoverCount = (crossoverIt != crossoverAllocation.crossoverCountPerSpecies.end()) ? 
                                 crossoverIt->second : 0;
        uint32_t guaranteedUnprotectedMutations = (guaranteedIt != equilibriumMutationAllocation.guaranteedUnprotectedMutationsPerSpecies.end()) ? 
                                                  guaranteedIt->second : 0;
        
        // Calculate remaining genomes after all prior allocations
        // For species below equilibrium, we allow expansion to equilibrium target
        uint32_t allocatedSoFar = eliteCount + crossoverCount + guaranteedUnprotectedMutations;
        uint32_t targetSize = std::max(data.instructionSetsSize, healthAssessment.equilibriumTarget);
        uint32_t remainingGenomes = (targetSize > allocatedSoFar) ? 
                                   (targetSize - allocatedSoFar) : 0;
        
        // Apply protection split to remaining genomes using exact integer arithmetic
        float protectionPercentage = params.getProtectionPercentage();
        
        // Calculate protected remainder using deterministic integer arithmetic
        // This ensures unprotected + protected = remaining with consistent behavior
        uint32_t protectedRemainder;
        if (protectionPercentage <= 0.0f) {
            protectedRemainder = 0;
        } else if (protectionPercentage >= 1.0f) {
            protectedRemainder = remainingGenomes;
        } else {
            // Use integer multiplication and division to avoid floating point precision issues
            // Formula: protectedRemainder = ceil(remainingGenomes * protectionPercentage)
            // Implemented as: (remainingGenomes * percentage_numerator + percentage_denominator - 1) / percentage_denominator
            
            // Convert percentage to rational representation with reasonable precision
            constexpr uint32_t PERCENTAGE_PRECISION = 10000; // 0.01% precision
            uint32_t protectionNumerator = static_cast<uint32_t>(protectionPercentage * PERCENTAGE_PRECISION + 0.5f);
            
            // Calculate using integer arithmetic with proper rounding
            uint64_t protectedProduct = static_cast<uint64_t>(remainingGenomes) * protectionNumerator;
            protectedRemainder = static_cast<uint32_t>((protectedProduct + PERCENTAGE_PRECISION - 1) / PERCENTAGE_PRECISION);
            
            // Clamp to valid range
            protectedRemainder = std::min(protectedRemainder, remainingGenomes);
        }
        
        uint32_t unprotectedRemainder = remainingGenomes - protectedRemainder;
        
        // Store results
        result.unprotectedRemainderPerSpecies[speciesId] = unprotectedRemainder;
        result.protectedRemainderPerSpecies[speciesId] = protectedRemainder;
        result.totalUnprotectedRemainder += unprotectedRemainder;
        result.totalProtectedRemainder += protectedRemainder;
        
        #ifdef DEBUG
        // Validation
        assert(allocatedSoFar <= targetSize);
        assert(unprotectedRemainder + protectedRemainder == remainingGenomes);
        assert(allocatedSoFar + remainingGenomes == targetSize);
        assert(eliteCount + crossoverCount + guaranteedUnprotectedMutations + 
               unprotectedRemainder + protectedRemainder == targetSize);
        #endif
    }
    
    return result;
}

// Phase 6: Instruction Set Generation with exact count matching
GenerationInstructionSets generateInstructionSets(
    const std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const EliteAllocationResult& eliteAllocation,
    const CrossoverAllocationResult& crossoverAllocation,
    const EquilibriumMutationAllocationResult& equilibriumMutationAllocation,
    const ProtectionAssignmentResult& protectionAssignment,
    const GenerationPlannerParams& params
) {
    auto logger = LOGGER("population.GenerationPlanner");
    GenerationInstructionSets instructionSets;
    
    for (const auto& [speciesId, data] : speciesData) {
        if (data.isMarkedForElimination || data.currentPopulationSize == 0) {
            // Eliminated species or empty species get empty instruction sets
            instructionSets[speciesId] = SpeciesInstructionSet{};
            continue;
        }
        
        // Get allocation counts for this species
        auto eliteIt = eliteAllocation.eliteCountPerSpecies.find(speciesId);
        auto crossoverIt = crossoverAllocation.crossoverCountPerSpecies.find(speciesId);
        auto guaranteedIt = equilibriumMutationAllocation.guaranteedUnprotectedMutationsPerSpecies.find(speciesId);
        auto unprotectedIt = protectionAssignment.unprotectedRemainderPerSpecies.find(speciesId);
        auto protectedIt = protectionAssignment.protectedRemainderPerSpecies.find(speciesId);
        
        uint32_t eliteCount = (eliteIt != eliteAllocation.eliteCountPerSpecies.end()) ? eliteIt->second : 0;
        uint32_t crossoverCount = (crossoverIt != crossoverAllocation.crossoverCountPerSpecies.end()) ? crossoverIt->second : 0;
        uint32_t guaranteedUnprotected = (guaranteedIt != equilibriumMutationAllocation.guaranteedUnprotectedMutationsPerSpecies.end()) ? guaranteedIt->second : 0;
        uint32_t unprotectedRemainder = (unprotectedIt != protectionAssignment.unprotectedRemainderPerSpecies.end()) ? unprotectedIt->second : 0;
        uint32_t protectedRemainder = (protectedIt != protectionAssignment.protectedRemainderPerSpecies.end()) ? protectedIt->second : 0;
        
        // Calculate total planned allocation for this species
        uint32_t totalPlannedAllocation = eliteCount + crossoverCount + guaranteedUnprotected + unprotectedRemainder + protectedRemainder;
        
        // Initialize instruction set for this species
        SpeciesInstructionSet& speciesInstructions = instructionSets[speciesId];
        speciesInstructions.reserve(totalPlannedAllocation);
        
        LOG_DEBUG("Species {}: Planning {} total instructions (elite={}, crossover={}, guaranteedUnprotected={}, unprotectedRemainder={}, protectedRemainder={})", 
                 speciesId, totalPlannedAllocation, eliteCount, crossoverCount, guaranteedUnprotected, unprotectedRemainder, protectedRemainder);
        
        uint32_t currentPriority = 0;
        
        // 1. Generate PRESERVE instructions (highest priority)
        for (uint32_t i = 0; i < eliteCount; ++i) {
            uint32_t relativeParentIndex = i; // Elite preserve top performers (relative indices 0, 1, 2...)
            speciesInstructions.push_back(
                ReproductiveInstruction::preserve(relativeParentIndex, currentPriority++)
            );
        }
        
        // 2. Generate CROSSOVER instructions 
        for (uint32_t i = 0; i < crossoverCount; ++i) {
            // Use proper fitness-proportional random selection from top half of species
            auto [parent1, parent2] = params.selectCrossoverParents(data.instructionSetsSize);
            
            speciesInstructions.push_back(
                ReproductiveInstruction::crossover(parent1, parent2, 0.7f, currentPriority++)
            );
        }
        
        // 3. Generate MUTATE_UNPROTECTED instructions (guaranteed + remainder)
        uint32_t totalUnprotectedMutations = guaranteedUnprotected + unprotectedRemainder;
        for (uint32_t i = 0; i < totalUnprotectedMutations; ++i) {
            // Fitness-proportional selection from entire species
            uint32_t relativeParentIndex = i % data.instructionSetsSize;
            speciesInstructions.push_back(
                ReproductiveInstruction::mutateUnprotected(relativeParentIndex, 0.1f, currentPriority++)
            );
        }
        
        // 4. Generate MUTATE_PROTECTED instructions (protected remainder)
        for (uint32_t i = 0; i < protectedRemainder; ++i) {
            // Protected mutations target worse performers (higher relative indices)
            uint32_t relativeParentIndex = (data.instructionSetsSize / 2) + (i % (data.instructionSetsSize - data.instructionSetsSize / 2));
            speciesInstructions.push_back(
                ReproductiveInstruction::mutateProtected(relativeParentIndex, 0.05f, currentPriority++)
            );
        }
        
        #ifdef DEBUG
        // Validate instruction count matches planned allocation
        assert(speciesInstructions.size() == totalPlannedAllocation);
        assert(speciesInstructions.size() == eliteCount + crossoverCount + totalUnprotectedMutations + protectedRemainder);
        
        // Validate all relative parent indices are within species bounds
        for (const auto& instruction : speciesInstructions) {
            assert(instruction.areParentIndicesValid(data.instructionSetsSize));
            assert(instruction.isValid());
        }
        #endif
    }
    
    return instructionSets;
}

// Main GenerationPlanner operator function
GenerationInstructionSets generationPlanner(
    std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
    const GenerationPlannerParams& params
) {
    auto logger = LOGGER("population.GenerationPlanner");
    LOG_DEBUG("generationPlanner: Starting with {} species", speciesData.size());
    
    // Phase 1: Population Health Assessment
    PopulationHealthAssessment healthAssessment = assessPopulationHealth(speciesData, params);
    
    // Phase 2: Elite Preservation Allocation
    EliteAllocationResult eliteAllocation = allocateElitePreservation(speciesData, params);
    
    // Phase 3: Crossover Allocation
    CrossoverAllocationResult crossoverAllocation = allocateCrossover(speciesData, eliteAllocation, params);
    
    // Phase 4: Equilibrium-Driven Mutation Allocation
    EquilibriumMutationAllocationResult equilibriumMutationAllocation = allocateEquilibriumMutations(
        speciesData, healthAssessment, eliteAllocation, crossoverAllocation, params);
    
    // Phase 5: Protection Assignment for Remaining Genomes
    ProtectionAssignmentResult protectionAssignment = assignProtectionForRemainingGenomes(
        speciesData, healthAssessment, eliteAllocation, crossoverAllocation, equilibriumMutationAllocation, params);
    
    // Phase 6: Instruction Set Generation
    GenerationInstructionSets instructionSets = generateInstructionSets(
        speciesData, eliteAllocation, crossoverAllocation, equilibriumMutationAllocation, protectionAssignment, params);
    
    // Log final instruction set totals
    uint32_t totalInstructionSets = 0;
    for (const auto& [speciesId, instructions] : instructionSets) {
        totalInstructionSets += instructions.size();
        LOG_DEBUG("Final species {}: {} instruction sets", speciesId, instructions.size());
    }
    LOG_DEBUG("Total instruction sets generated: {}", totalInstructionSets);
    
    #ifdef DEBUG
    // Final validation: verify atomic instruction generation guarantee
    for (const auto& [speciesId, data] : speciesData) {
        auto instructionIt = instructionSets.find(speciesId);
        assert(instructionIt != instructionSets.end());
        
        if (data.isMarkedForElimination || data.instructionSetsSize == 0) {
            assert(instructionIt->second.empty());
        } else {
            assert(instructionIt->second.size() == data.instructionSetsSize);
        }
    }
    
    // Validate total instruction count matches total instruction set sizes
    uint32_t totalInstructionsGenerated = 0;
    uint32_t totalInstructionSetsExpected = 0;
    for (const auto& [speciesId, data] : speciesData) {
        if (!data.isMarkedForElimination) {
            totalInstructionSetsExpected += data.instructionSetsSize;
            totalInstructionsGenerated += instructionSets[speciesId].size();
        }
    }
    assert(totalInstructionsGenerated == totalInstructionSetsExpected);
    #endif
    
    return instructionSets;
}

} // namespace Population