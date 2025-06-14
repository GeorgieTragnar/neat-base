#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

namespace Population {

// Enumeration for reproductive operation types
enum class OperationType {
    PRESERVE,           // Direct copy without evolution (elite preservation)
    CROSSOVER,          // Genetic recombination between two parents
    MUTATE_UNPROTECTED, // Evolutionary mutation without protection tracking
    MUTATE_PROTECTED    // Conservative mutation with protection tracking
};

// Evolution parameters for mutation and crossover operations
struct EvolutionParameters {
    float mutationRate = 0.1f;           // Base mutation rate for MUTATE operations
    float crossoverRate = 0.7f;          // Crossover probability for CROSSOVER operations
    bool enableProtectionTracking = false; // Whether to track protection for this instruction
    uint32_t priority = 0;               // Execution ordering hint (lower = higher priority)
};

// Reproductive instruction specifying how to generate a new genome
struct ReproductiveInstruction {
    OperationType operationType;         // Type of reproductive operation to perform
    std::vector<uint32_t> relativeParentIndices; // Parent indices relative to species (0-based)
    std::vector<uint32_t> globalParentIndices; // Global parent indices into the master multimap with fitness results
    EvolutionParameters evolutionParams; // Parameters for evolution operations
    
    // Constructor for PRESERVE operations (single parent, no evolution)
    static ReproductiveInstruction preserve(uint32_t relativeParentIndex, uint32_t priority = 0) {
        return ReproductiveInstruction{
            OperationType::PRESERVE,
            {relativeParentIndex},
            {},
            {0.0f, 0.0f, false, priority}
        };
    }
    
    // Constructor for CROSSOVER operations (two parents)
    static ReproductiveInstruction crossover(
        uint32_t relativeParent1Index, 
        uint32_t relativeParent2Index,
        float crossoverRate = 0.7f,
        uint32_t priority = 0
    ) {
        return ReproductiveInstruction{
            OperationType::CROSSOVER,
            {relativeParent1Index, relativeParent2Index},
            {},
            {0.0f, crossoverRate, false, priority}
        };
    }
    
    // Constructor for MUTATE_UNPROTECTED operations (single parent with mutation)
    static ReproductiveInstruction mutateUnprotected(
        uint32_t relativeParentIndex,
        float mutationRate = 0.1f,
        uint32_t priority = 0
    ) {
        return ReproductiveInstruction{
            OperationType::MUTATE_UNPROTECTED,
            {relativeParentIndex},
            {},
            {mutationRate, 0.0f, false, priority}
        };
    }
    
    // Constructor for MUTATE_PROTECTED operations (single parent with protection tracking)
    static ReproductiveInstruction mutateProtected(
        uint32_t relativeParentIndex,
        float mutationRate = 0.05f, // Lower default mutation rate for protected
        uint32_t priority = 0
    ) {
        return ReproductiveInstruction{
            OperationType::MUTATE_PROTECTED,
            {relativeParentIndex},
            {},
            {mutationRate, 0.0f, true, priority}
        };
    }
    
    // Validation methods
    bool isValid() const {
        switch (operationType) {
            case OperationType::PRESERVE:
            case OperationType::MUTATE_UNPROTECTED:
            case OperationType::MUTATE_PROTECTED:
                return relativeParentIndices.size() == 1;
            case OperationType::CROSSOVER:
                return relativeParentIndices.size() == 2 && 
                       relativeParentIndices[0] != relativeParentIndices[1];
            default:
                return false;
        }
    }
    
    // Check if this instruction requires protection tracking
    bool requiresProtectionTracking() const {
        return evolutionParams.enableProtectionTracking;
    }
    
    // Get the number of parents required for this operation
    size_t getParentCount() const {
        return relativeParentIndices.size();
    }
    
    // Check if all parent indices are within valid range for species size
    bool areParentIndicesValid(uint32_t speciesSize) const {
        if (speciesSize == 0) return false;
        
        for (uint32_t index : relativeParentIndices) {
            if (index >= speciesSize) {
                return false;
            }
        }
        return true;
    }
};

// Collection of reproductive instructions for a species
using SpeciesInstructionSet = std::vector<ReproductiveInstruction>;

// Complete instruction sets for all species
using GenerationInstructionSets = std::unordered_map<uint32_t, SpeciesInstructionSet>;

} // namespace Population