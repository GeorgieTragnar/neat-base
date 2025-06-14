#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <random>

namespace Population {

// Behavior options for rank scaling when species rank exceeds array bounds
enum class ScalingBoundaryBehavior {
    USE_MINIMUM,        // Use minimum value from array (current behavior)
    USE_LAST_VALUE,     // Repeat the last value in the array
    LINEAR_DECAY,       // Continue with linear decay from last value
    EXPONENTIAL_DECAY,  // Continue with exponential decay from last value
    CONSTANT_MINIMUM    // Use a constant minimum value
};

// Forward declarations
struct DynamicSpeciesData;
struct ReproductiveInstruction;
using SpeciesInstructionSet = std::vector<ReproductiveInstruction>;
using GenerationInstructionSets = std::unordered_map<uint32_t, SpeciesInstructionSet>;

// Configuration parameters for the GenerationPlanner operator
class GenerationPlannerParams {
public:
    // Delete default constructor - all parameters must be explicitly set
    GenerationPlannerParams() = delete;
    
    // Constructor with all required parameters
    GenerationPlannerParams(
        // Elite Preservation Parameters
        uint32_t baseEliteCount,
        const std::vector<float>& eliteScalingByRank,
        float maxElitePercentage,
        ScalingBoundaryBehavior eliteBoundaryBehavior,
        float eliteConstantMinimum,
        float eliteDecayRate,
        
        // Crossover Allocation Parameters
        uint32_t baseCrossoverSlots,
        const std::vector<float>& crossoverScalingByRank,
        uint32_t minSpeciesSizeForCrossover,
        ScalingBoundaryBehavior crossoverBoundaryBehavior,
        float crossoverConstantMinimum,
        float crossoverDecayRate,
        
        // Equilibrium Parameters
        uint32_t targetTotalPopulation,
        float equilibriumBias,
        uint32_t minSpeciesSize,
        
        // Protection Parameters
        float protectionPercentage,
        uint32_t protectionThreshold,
        uint32_t speciesEliminationThreshold,
        
        // Random Number Generation
        uint32_t randomSeed
    );

protected:
    // Friend declaration for the main operator function
    friend GenerationInstructionSets generationPlanner(
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const GenerationPlannerParams& params
    );

    // Elite Preservation Parameters
    const uint32_t _baseEliteCount;           // Minimum elite preservation per species
    const std::vector<float> _eliteScalingByRank; // Rank-based multipliers [rank1: 3.0x, rank2: 2.5x, ...]
    const float _maxElitePercentage;          // Cap on elite percentage per species
    const ScalingBoundaryBehavior _eliteBoundaryBehavior; // How to handle ranks beyond array bounds
    const float _eliteConstantMinimum;        // Constant minimum for CONSTANT_MINIMUM behavior
    const float _eliteDecayRate;              // Decay rate for decay behaviors

    // Crossover Allocation Parameters  
    const uint32_t _baseCrossoverSlots;       // Base crossover allocation per species
    const std::vector<float> _crossoverScalingByRank; // Rank-based crossover multipliers
    const uint32_t _minSpeciesSizeForCrossover; // Minimum instructionSetSize to enable crossover
    const ScalingBoundaryBehavior _crossoverBoundaryBehavior; // How to handle ranks beyond array bounds
    const float _crossoverConstantMinimum;    // Constant minimum for CONSTANT_MINIMUM behavior
    const float _crossoverDecayRate;          // Decay rate for decay behaviors

    // Equilibrium Parameters
    const uint32_t _targetTotalPopulation;    // Desired total instruction sets across all species
    const float _equilibriumBias;             // Strength of equilibrium force (0.0 = no bias, 1.0 = maximum)
    const uint32_t _minSpeciesSize;           // Absolute minimum instruction sets per species

    // Protection Parameters
    const float _protectionPercentage;        // Percentage of remaining genomes marked for protection
    const uint32_t _protectionThreshold;      // Protection counter limit before instruction set denial
    const uint32_t _speciesEliminationThreshold; // Species protection rating limit

    // Random Number Generation (mutable for state-changing operations)
    mutable std::mt19937 _randomGenerator;    // Mersenne Twister generator for parent selection

public:
    // Accessor methods for parameter validation and debugging
    
    // Elite parameters
    uint32_t getBaseEliteCount() const { return _baseEliteCount; }
    const std::vector<float>& getEliteScalingByRank() const { return _eliteScalingByRank; }
    float getMaxElitePercentage() const { return _maxElitePercentage; }
    
    // Crossover parameters
    uint32_t getBaseCrossoverSlots() const { return _baseCrossoverSlots; }
    const std::vector<float>& getCrossoverScalingByRank() const { return _crossoverScalingByRank; }
    uint32_t getMinSpeciesSizeForCrossover() const { return _minSpeciesSizeForCrossover; }
    
    // Equilibrium parameters
    uint32_t getTargetTotalPopulation() const { return _targetTotalPopulation; }
    float getEquilibriumBias() const { return _equilibriumBias; }
    uint32_t getMinSpeciesSize() const { return _minSpeciesSize; }
    
    // Protection parameters
    float getProtectionPercentage() const { return _protectionPercentage; }
    uint32_t getProtectionThreshold() const { return _protectionThreshold; }
    uint32_t getSpeciesEliminationThreshold() const { return _speciesEliminationThreshold; }
    
    // Utility methods for scaling factor retrieval
    float getEliteScalingFactor(uint32_t rank) const;
    float getCrossoverScalingFactor(uint32_t rank) const;
    
    // Validation methods
    bool isValidRank(uint32_t rank) const;
    bool isValidSpeciesSize(uint32_t size) const;
    bool canEnableCrossover(uint32_t speciesSize) const;
    
    // Random number generation for parent selection
    std::mt19937& getRandomGenerator() const;
    
    // Helper method for fitness-proportional parent selection from top half of species
    std::pair<uint32_t, uint32_t> selectCrossoverParents(uint32_t speciesSize) const;

private:
    // Helper method to calculate scaling factors for ranks beyond array bounds
    float calculateBoundaryScalingFactor(const std::vector<float>& scalingArray, 
                                       ScalingBoundaryBehavior behavior,
                                       float constantMinimum, float decayRate, 
                                       uint32_t rank) const;
};

// Factory function for creating common parameter configurations
class GenerationPlannerParamsFactory {
public:
    // Conservative configuration - emphasizes stability and diversity preservation
    static GenerationPlannerParams createConservative(uint32_t targetPopulation, uint32_t randomSeed);
    
    // Aggressive configuration - emphasizes performance and rapid convergence
    static GenerationPlannerParams createAggressive(uint32_t targetPopulation, uint32_t randomSeed);
    
    // Balanced configuration - balanced exploration and exploitation
    static GenerationPlannerParams createBalanced(uint32_t targetPopulation, uint32_t randomSeed);
};

} // namespace Population