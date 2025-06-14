#include "GenerationPlannerParams.hpp"

namespace Population {

// GenerationPlannerParams constructor
GenerationPlannerParams::GenerationPlannerParams(
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
) : _baseEliteCount(baseEliteCount),
    _eliteScalingByRank(eliteScalingByRank),
    _maxElitePercentage(maxElitePercentage),
    _eliteBoundaryBehavior(eliteBoundaryBehavior),
    _eliteConstantMinimum(eliteConstantMinimum),
    _eliteDecayRate(eliteDecayRate),
    _baseCrossoverSlots(baseCrossoverSlots),
    _crossoverScalingByRank(crossoverScalingByRank),
    _minSpeciesSizeForCrossover(minSpeciesSizeForCrossover),
    _crossoverBoundaryBehavior(crossoverBoundaryBehavior),
    _crossoverConstantMinimum(crossoverConstantMinimum),
    _crossoverDecayRate(crossoverDecayRate),
    _targetTotalPopulation(targetTotalPopulation),
    _equilibriumBias(equilibriumBias),
    _minSpeciesSize(minSpeciesSize),
    _protectionPercentage(protectionPercentage),
    _protectionThreshold(protectionThreshold),
    _speciesEliminationThreshold(speciesEliminationThreshold),
    _randomGenerator(randomSeed) {
    
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
    
    // Validate boundary behavior parameters
    assert(eliteConstantMinimum >= 0.0f && "Elite constant minimum must be non-negative");
    assert(eliteDecayRate > 0.0f && eliteDecayRate <= 1.0f && "Elite decay rate must be in (0.0, 1.0]");
    assert(crossoverConstantMinimum >= 0.0f && "Crossover constant minimum must be non-negative");
    assert(crossoverDecayRate > 0.0f && crossoverDecayRate <= 1.0f && "Crossover decay rate must be in (0.0, 1.0]");
}

// Utility methods for scaling factor retrieval
float GenerationPlannerParams::getEliteScalingFactor(uint32_t rank) const {
    if (rank == 0) {
        // Invalid rank (ranks are 1-based)
        return 0.0f;
    }
    if (rank <= _eliteScalingByRank.size()) {
        // Rank within array bounds
        return _eliteScalingByRank[rank - 1]; // Convert 1-based rank to 0-based index
    }
    // Rank beyond array bounds - apply boundary behavior
    return calculateBoundaryScalingFactor(_eliteScalingByRank, _eliteBoundaryBehavior, 
                                        _eliteConstantMinimum, _eliteDecayRate, rank);
}

float GenerationPlannerParams::getCrossoverScalingFactor(uint32_t rank) const {
    if (rank == 0) {
        // Invalid rank (ranks are 1-based)
        return 0.0f;
    }
    if (rank <= _crossoverScalingByRank.size()) {
        // Rank within array bounds
        return _crossoverScalingByRank[rank - 1]; // Convert 1-based rank to 0-based index
    }
    // Rank beyond array bounds - apply boundary behavior
    return calculateBoundaryScalingFactor(_crossoverScalingByRank, _crossoverBoundaryBehavior, 
                                        _crossoverConstantMinimum, _crossoverDecayRate, rank);
}

// Validation methods
bool GenerationPlannerParams::isValidRank(uint32_t rank) const {
    return rank > 0; // Ranks are 1-based
}

bool GenerationPlannerParams::isValidSpeciesSize(uint32_t size) const {
    return size >= _minSpeciesSize;
}

bool GenerationPlannerParams::canEnableCrossover(uint32_t speciesSize) const {
    return speciesSize >= _minSpeciesSizeForCrossover;
}

// Random number generation for parent selection
std::mt19937& GenerationPlannerParams::getRandomGenerator() const {
    return _randomGenerator;
}

// Helper method for fitness-proportional parent selection from top half of species
std::pair<uint32_t, uint32_t> GenerationPlannerParams::selectCrossoverParents(uint32_t speciesSize) const {
    if (speciesSize < 2) {
        // Degenerate case - return same parent twice (shouldn't happen in practice)
        return {0, 0};
    }
    
    // Select from top half of species (better performers have lower indices)
    uint32_t maxParentIndex = std::max(2u, speciesSize / 2);
    
    // Use weighted selection favoring better performers (lower indices)
    // Weight = maxParentIndex - index (so index 0 has highest weight)
    std::vector<uint32_t> weights;
    weights.reserve(maxParentIndex);
    for (uint32_t i = 0; i < maxParentIndex; ++i) {
        weights.push_back(maxParentIndex - i);
    }
    
    // Create distribution for weighted selection
    std::discrete_distribution<uint32_t> distribution(weights.begin(), weights.end());
    
    // Select two different parents
    uint32_t parent1 = distribution(_randomGenerator);
    uint32_t parent2;
    
    // Ensure different parents
    do {
        parent2 = distribution(_randomGenerator);
    } while (parent2 == parent1 && maxParentIndex > 1);
    
    return {parent1, parent2};
}

// Helper method to calculate scaling factors for ranks beyond array bounds
float GenerationPlannerParams::calculateBoundaryScalingFactor(const std::vector<float>& scalingArray, 
                                   ScalingBoundaryBehavior behavior,
                                   float constantMinimum, float decayRate, 
                                   uint32_t rank) const {
    switch (behavior) {
        case ScalingBoundaryBehavior::USE_MINIMUM:
            return *std::min_element(scalingArray.begin(), scalingArray.end());
            
        case ScalingBoundaryBehavior::USE_LAST_VALUE:
            return scalingArray.back();
            
        case ScalingBoundaryBehavior::LINEAR_DECAY: {
            if (scalingArray.size() < 2) {
                return scalingArray.back();
            }
            // Calculate linear decay rate from last two values
            float lastValue = scalingArray.back();
            float secondLastValue = scalingArray[scalingArray.size() - 2];
            float linearDecayStep = secondLastValue - lastValue;
            
            // Apply decay for ranks beyond array
            uint32_t excessRanks = rank - static_cast<uint32_t>(scalingArray.size());
            float decayedValue = lastValue - (linearDecayStep * excessRanks);
            
            // Ensure non-negative result
            return std::max(0.0f, decayedValue);
        }
        
        case ScalingBoundaryBehavior::EXPONENTIAL_DECAY: {
            float lastValue = scalingArray.back();
            uint32_t excessRanks = rank - static_cast<uint32_t>(scalingArray.size());
            
            // Apply exponential decay: value = lastValue * (decayRate ^ excessRanks)
            float decayedValue = lastValue;
            for (uint32_t i = 0; i < excessRanks; ++i) {
                decayedValue *= decayRate;
            }
            
            return decayedValue;
        }
        
        case ScalingBoundaryBehavior::CONSTANT_MINIMUM:
            return constantMinimum;
            
        default:
            // Fallback to minimum behavior
            return *std::min_element(scalingArray.begin(), scalingArray.end());
    }
}

// Factory method implementations
GenerationPlannerParams GenerationPlannerParamsFactory::createConservative(uint32_t targetPopulation, uint32_t randomSeed) {
    return GenerationPlannerParams(
        // Elite parameters - conservative preservation
        2,                                    // baseEliteCount
        {2.0f, 1.8f, 1.6f, 1.4f, 1.2f},     // eliteScalingByRank (moderate scaling)
        0.5f,                                 // maxElitePercentage (up to 50% can be elite)
        ScalingBoundaryBehavior::LINEAR_DECAY, // eliteBoundaryBehavior (gradual decay)
        0.1f,                                 // eliteConstantMinimum
        0.9f,                                 // eliteDecayRate
        
        // Crossover parameters - moderate crossover
        1,                                    // baseCrossoverSlots
        {1.5f, 1.3f, 1.1f, 1.0f, 0.8f},     // crossoverScalingByRank
        3,                                    // minSpeciesSizeForCrossover
        ScalingBoundaryBehavior::USE_LAST_VALUE, // crossoverBoundaryBehavior (stable)
        0.0f,                                 // crossoverConstantMinimum
        0.8f,                                 // crossoverDecayRate
        
        // Equilibrium parameters - strong equilibrium force
        targetPopulation,                     // targetTotalPopulation
        0.8f,                                 // equilibriumBias (strong equilibrium)
        2,                                    // minSpeciesSize
        
        // Protection parameters - high protection
        0.4f,                                 // protectionPercentage (40% protected)
        3,                                    // protectionThreshold
        5,                                    // speciesEliminationThreshold
        
        // Random number generation
        randomSeed                            // randomSeed
    );
}

GenerationPlannerParams GenerationPlannerParamsFactory::createAggressive(uint32_t targetPopulation, uint32_t randomSeed) {
    return GenerationPlannerParams(
        // Elite parameters - aggressive preservation of top performers
        1,                                    // baseEliteCount
        {4.0f, 3.0f, 2.0f, 1.5f, 1.0f},     // eliteScalingByRank (strong scaling)
        0.3f,                                 // maxElitePercentage (up to 30% can be elite)
        ScalingBoundaryBehavior::EXPONENTIAL_DECAY, // eliteBoundaryBehavior (rapid decay)
        0.05f,                                // eliteConstantMinimum
        0.7f,                                 // eliteDecayRate
        
        // Crossover parameters - high crossover activity
        2,                                    // baseCrossoverSlots
        {3.0f, 2.5f, 2.0f, 1.5f, 1.0f},     // crossoverScalingByRank
        2,                                    // minSpeciesSizeForCrossover
        ScalingBoundaryBehavior::EXPONENTIAL_DECAY, // crossoverBoundaryBehavior (rapid decay)
        0.0f,                                 // crossoverConstantMinimum
        0.6f,                                 // crossoverDecayRate
        
        // Equilibrium parameters - moderate equilibrium force
        targetPopulation,                     // targetTotalPopulation
        0.5f,                                 // equilibriumBias (moderate equilibrium)
        1,                                    // minSpeciesSize
        
        // Protection parameters - low protection (high selection pressure)
        0.2f,                                 // protectionPercentage (20% protected)
        2,                                    // protectionThreshold
        3,                                    // speciesEliminationThreshold
        
        // Random number generation
        randomSeed                            // randomSeed
    );
}

GenerationPlannerParams GenerationPlannerParamsFactory::createBalanced(uint32_t targetPopulation, uint32_t randomSeed) {
    return GenerationPlannerParams(
        // Elite parameters - balanced preservation
        1,                                    // baseEliteCount
        {2.5f, 2.0f, 1.7f, 1.4f, 1.1f},     // eliteScalingByRank
        0.4f,                                 // maxElitePercentage (up to 40% can be elite)
        ScalingBoundaryBehavior::LINEAR_DECAY, // eliteBoundaryBehavior (balanced decay)
        0.1f,                                 // eliteConstantMinimum
        0.8f,                                 // eliteDecayRate
        
        // Crossover parameters - balanced crossover
        1,                                    // baseCrossoverSlots
        {2.0f, 1.7f, 1.4f, 1.2f, 1.0f},     // crossoverScalingByRank
        2,                                    // minSpeciesSizeForCrossover
        ScalingBoundaryBehavior::LINEAR_DECAY, // crossoverBoundaryBehavior (balanced decay)
        0.1f,                                 // crossoverConstantMinimum
        0.8f,                                 // crossoverDecayRate
        
        // Equilibrium parameters - balanced equilibrium force
        targetPopulation,                     // targetTotalPopulation
        0.6f,                                 // equilibriumBias (moderate-strong equilibrium)
        1,                                    // minSpeciesSize
        
        // Protection parameters - moderate protection
        0.3f,                                 // protectionPercentage (30% protected)
        3,                                    // protectionThreshold
        4,                                    // speciesEliminationThreshold
        
        // Random number generation
        randomSeed                            // randomSeed
    );
}

} // namespace Population