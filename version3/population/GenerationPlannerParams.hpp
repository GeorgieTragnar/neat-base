#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <cstdint>

namespace Population {

// Forward declaration
struct DynamicSpeciesData;

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
        
        // Crossover Allocation Parameters
        uint32_t baseCrossoverSlots,
        const std::vector<float>& crossoverScalingByRank,
        uint32_t minSpeciesSizeForCrossover,
        
        // Equilibrium Parameters
        uint32_t targetTotalPopulation,
        float equilibriumBias,
        uint32_t minSpeciesSize,
        
        // Protection Parameters
        float protectionPercentage,
        uint32_t protectionThreshold,
        uint32_t speciesEliminationThreshold
    ) : _baseEliteCount(baseEliteCount),
        _eliteScalingByRank(eliteScalingByRank),
        _maxElitePercentage(maxElitePercentage),
        _baseCrossoverSlots(baseCrossoverSlots),
        _crossoverScalingByRank(crossoverScalingByRank),
        _minSpeciesSizeForCrossover(minSpeciesSizeForCrossover),
        _targetTotalPopulation(targetTotalPopulation),
        _equilibriumBias(equilibriumBias),
        _minSpeciesSize(minSpeciesSize),
        _protectionPercentage(protectionPercentage),
        _protectionThreshold(protectionThreshold),
        _speciesEliminationThreshold(speciesEliminationThreshold) {
        
        // Parameter validation
        assert(baseEliteCount > 0 && "Base elite count must be positive");
        assert(!eliteScalingByRank.empty() && "Elite scaling array cannot be empty");
        assert(maxElitePercentage >= 0.0f && maxElitePercentage <= 1.0f && "Max elite percentage must be [0.0, 1.0]");
        
        assert(baseCrossoverSlots >= 0 && "Base crossover slots must be non-negative");
        assert(!crossoverScalingByRank.empty() && "Crossover scaling array cannot be empty");
        assert(minSpeciesSizeForCrossover >= 2 && "Min species size for crossover must be at least 2");
        
        assert(targetTotalPopulation > 0 && "Target total population must be positive");
        assert(equilibriumBias >= 0.0f && equilibriumBias <= 1.0f && "Equilibrium bias must be [0.0, 1.0]");
        assert(minSpeciesSize > 0 && "Min species size must be positive");
        
        assert(protectionPercentage >= 0.0f && protectionPercentage <= 1.0f && "Protection percentage must be [0.0, 1.0]");
        assert(protectionThreshold > 0 && "Protection threshold must be positive");
        assert(speciesEliminationThreshold > 0 && "Species elimination threshold must be positive");
        
        // Validate scaling arrays have positive values
        for (size_t i = 0; i < eliteScalingByRank.size(); ++i) {
            assert(eliteScalingByRank[i] > 0.0f && "Elite scaling factors must be positive");
        }
        for (size_t i = 0; i < crossoverScalingByRank.size(); ++i) {
            assert(crossoverScalingByRank[i] >= 0.0f && "Crossover scaling factors must be non-negative");
        }
    }

protected:
    // Friend declaration for the main operator function
    friend void generationPlanner(
        std::unordered_map<uint32_t, DynamicSpeciesData>& speciesData,
        const GenerationPlannerParams& params
    );

    // Elite Preservation Parameters
    const uint32_t _baseEliteCount;           // Minimum elite preservation per species
    const std::vector<float> _eliteScalingByRank; // Rank-based multipliers [rank1: 3.0x, rank2: 2.5x, ...]
    const float _maxElitePercentage;          // Cap on elite percentage per species

    // Crossover Allocation Parameters  
    const uint32_t _baseCrossoverSlots;       // Base crossover allocation per species
    const std::vector<float> _crossoverScalingByRank; // Rank-based crossover multipliers
    const uint32_t _minSpeciesSizeForCrossover; // Minimum instructionSetSize to enable crossover

    // Equilibrium Parameters
    const uint32_t _targetTotalPopulation;    // Desired total instruction sets across all species
    const float _equilibriumBias;             // Strength of equilibrium force (0.0 = no bias, 1.0 = maximum)
    const uint32_t _minSpeciesSize;           // Absolute minimum instruction sets per species

    // Protection Parameters
    const float _protectionPercentage;        // Percentage of remaining genomes marked for protection
    const uint32_t _protectionThreshold;      // Protection counter limit before instruction set denial
    const uint32_t _speciesEliminationThreshold; // Species protection rating limit

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
    float getEliteScalingFactor(uint32_t rank) const {
        if (rank == 0 || rank > _eliteScalingByRank.size()) {
            // Use minimum scaling factor for ranks beyond array or invalid rank
            return *std::min_element(_eliteScalingByRank.begin(), _eliteScalingByRank.end());
        }
        return _eliteScalingByRank[rank - 1]; // Convert 1-based rank to 0-based index
    }
    
    float getCrossoverScalingFactor(uint32_t rank) const {
        if (rank == 0 || rank > _crossoverScalingByRank.size()) {
            // Use minimum scaling factor for ranks beyond array or invalid rank
            return *std::min_element(_crossoverScalingByRank.begin(), _crossoverScalingByRank.end());
        }
        return _crossoverScalingByRank[rank - 1]; // Convert 1-based rank to 0-based index
    }
    
    // Validation methods
    bool isValidRank(uint32_t rank) const {
        return rank > 0; // Ranks are 1-based
    }
    
    bool isValidSpeciesSize(uint32_t size) const {
        return size >= _minSpeciesSize;
    }
    
    bool canEnableCrossover(uint32_t speciesSize) const {
        return speciesSize >= _minSpeciesSizeForCrossover;
    }
};

// Factory function for creating common parameter configurations
class GenerationPlannerParamsFactory {
public:
    // Conservative configuration - emphasizes stability and diversity preservation
    static GenerationPlannerParams createConservative(uint32_t targetPopulation) {
        return GenerationPlannerParams(
            // Elite parameters - conservative preservation
            2,                                    // baseEliteCount
            {2.0f, 1.8f, 1.6f, 1.4f, 1.2f},     // eliteScalingByRank (moderate scaling)
            0.5f,                                 // maxElitePercentage (up to 50% can be elite)
            
            // Crossover parameters - moderate crossover
            1,                                    // baseCrossoverSlots
            {1.5f, 1.3f, 1.1f, 1.0f, 0.8f},     // crossoverScalingByRank
            3,                                    // minSpeciesSizeForCrossover
            
            // Equilibrium parameters - strong equilibrium force
            targetPopulation,                     // targetTotalPopulation
            0.8f,                                 // equilibriumBias (strong equilibrium)
            2,                                    // minSpeciesSize
            
            // Protection parameters - high protection
            0.4f,                                 // protectionPercentage (40% protected)
            3,                                    // protectionThreshold
            5                                     // speciesEliminationThreshold
        );
    }
    
    // Aggressive configuration - emphasizes performance and rapid convergence
    static GenerationPlannerParams createAggressive(uint32_t targetPopulation) {
        return GenerationPlannerParams(
            // Elite parameters - aggressive preservation of top performers
            1,                                    // baseEliteCount
            {4.0f, 3.0f, 2.0f, 1.5f, 1.0f},     // eliteScalingByRank (strong scaling)
            0.3f,                                 // maxElitePercentage (up to 30% can be elite)
            
            // Crossover parameters - high crossover activity
            2,                                    // baseCrossoverSlots
            {3.0f, 2.5f, 2.0f, 1.5f, 1.0f},     // crossoverScalingByRank
            2,                                    // minSpeciesSizeForCrossover
            
            // Equilibrium parameters - moderate equilibrium force
            targetPopulation,                     // targetTotalPopulation
            0.5f,                                 // equilibriumBias (moderate equilibrium)
            1,                                    // minSpeciesSize
            
            // Protection parameters - low protection (high selection pressure)
            0.2f,                                 // protectionPercentage (20% protected)
            2,                                    // protectionThreshold
            3                                     // speciesEliminationThreshold
        );
    }
    
    // Balanced configuration - balanced exploration and exploitation
    static GenerationPlannerParams createBalanced(uint32_t targetPopulation) {
        return GenerationPlannerParams(
            // Elite parameters - balanced preservation
            1,                                    // baseEliteCount
            {2.5f, 2.0f, 1.7f, 1.4f, 1.1f},     // eliteScalingByRank
            0.4f,                                 // maxElitePercentage (up to 40% can be elite)
            
            // Crossover parameters - balanced crossover
            1,                                    // baseCrossoverSlots
            {2.0f, 1.7f, 1.4f, 1.2f, 1.0f},     // crossoverScalingByRank
            2,                                    // minSpeciesSizeForCrossover
            
            // Equilibrium parameters - balanced equilibrium force
            targetPopulation,                     // targetTotalPopulation
            0.6f,                                 // equilibriumBias (moderate-strong equilibrium)
            1,                                    // minSpeciesSize
            
            // Protection parameters - moderate protection
            0.3f,                                 // protectionPercentage (30% protected)
            3,                                    // protectionThreshold
            4                                     // speciesEliminationThreshold
        );
    }
};

} // namespace Population