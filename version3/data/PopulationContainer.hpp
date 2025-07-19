#pragma once

#include <array>
#include <vector>
#include <map>
#include <cstdint>
#include <cassert>

#include "GlobalIndexRegistry.hpp"

#include "data_forward_declarations.inc"
#include "operator_forward_declarations.inc"

// Forward declarations
class Genome;

// Forward declarations
struct DynamicGenomeData;

template<typename FitnessResultType>
class PopulationContainer {
public:
    explicit PopulationContainer(GlobalIndexRegistry& registry) : _registry(registry) {}
    ~PopulationContainer() = default;

    // Public const accessors only
    const std::vector<Genome>& getGenomes(uint32_t generation) const {
        return _genomes[generation % 3];
    }
    
    const std::vector<DynamicGenomeData>& getGenomeData(uint32_t generation) const {
        return _genomeData[generation % 3];
    }
    
    const std::multimap<FitnessResultType, size_t>& getFitnessResults(uint32_t generation) const {
        return _fitnessResults[generation % 3];
    }
    
    // Const convenience methods for common access patterns with underflow protection
    const std::vector<Genome>& getCurrentGenomes(uint32_t currentGen) const {
        return getGenomes(currentGen);
    }
    
    const std::vector<Genome>& getLastGenomes(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation: currentGen underflow");
        return getGenomes(currentGen - 1);
    }
    
    const std::vector<Genome>& getGenerationBeforeLastGenomes(uint32_t currentGen) const {
        assert(currentGen > 1 && "Cannot access generation before last: currentGen underflow");
        return getGenomes(currentGen - 2);
    }
    
    const std::vector<DynamicGenomeData>& getCurrentGenomeData(uint32_t currentGen) const {
        return getGenomeData(currentGen);
    }
    
    const std::vector<DynamicGenomeData>& getLastGenomeData(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation data: currentGen underflow");
        return getGenomeData(currentGen - 1);
    }
    
    const std::vector<DynamicGenomeData>& getGenerationBeforeLastGenomeData(uint32_t currentGen) const {
        assert(currentGen > 1 && "Cannot access generation before last data: currentGen underflow");
        return getGenomeData(currentGen - 2);
    }
    
    const std::multimap<FitnessResultType, size_t>& getCurrentFitnessResults(uint32_t currentGen) const {
        return getFitnessResults(currentGen);
    }
    
    const std::multimap<FitnessResultType, size_t>& getLastFitnessResults(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation fitness: currentGen underflow");
        return getFitnessResults(currentGen - 1);
    }
    
    const std::multimap<FitnessResultType, size_t>& getGenerationBeforeLastFitnessResults(uint32_t currentGen) const {
        assert(currentGen > 1 && "Cannot access generation before last fitness: currentGen underflow");
        return getFitnessResults(currentGen - 2);
    }
    
    // Container operations - only fitness results get cleared, genomes/data handled via assignment
    void clearGenerationFitnessResults(uint32_t generation) {
        _fitnessResults[generation % 3].clear();
    }
    
    void reserveCapacity(size_t capacity) {
        for (size_t i = 0; i < 3; ++i) {
            _genomes[i].reserve(capacity);
            _genomeData[i].reserve(capacity);
        }
    }
    
    // Validation - check that all arrays within each generation have consistent sizes
    bool validateConsistency() const {
        for (size_t i = 0; i < 3; ++i) {
            if (_genomes[i].size() != _genomeData[i].size()) {
                return false;
            }
        }
        return true;
    }
    
    // Debugging/inspection methods
    size_t getGenerationSize(uint32_t generation) const {
        return _genomes[generation % 3].size();
    }
    
    size_t getGenerationCapacity(uint32_t generation) const {
        return _genomes[generation % 3].capacity();
    }

// protected:
    // Keep existing manual friend declaration
    friend class GlobalIndexRegistry;
    
#include "operator_template_friend_declarations.inc"

    // Mutable accessors for operators only
    std::vector<Genome>& getGenomes(uint32_t generation) {
        return _genomes[generation % 3];
    }
    
    std::vector<DynamicGenomeData>& getGenomeData(uint32_t generation) {
        return _genomeData[generation % 3];
    }
    
    std::multimap<FitnessResultType, size_t>& getFitnessResults(uint32_t generation) {
        return _fitnessResults[generation % 3];
    }
    
    // Mutable convenience methods for operators
    std::vector<Genome>& getCurrentGenomes(uint32_t currentGen) {
        return getGenomes(currentGen);
    }
    
    std::vector<Genome>& getLastGenomes(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation: currentGen underflow");
        return getGenomes(currentGen - 1);
    }
    
    std::vector<Genome>& getGenerationBeforeLastGenomes(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last: currentGen underflow");
        return getGenomes(currentGen - 2);
    }
    
    std::vector<DynamicGenomeData>& getCurrentGenomeData(uint32_t currentGen) {
        return getGenomeData(currentGen);
    }
    
    std::vector<DynamicGenomeData>& getLastGenomeData(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation data: currentGen underflow");
        return getGenomeData(currentGen - 1);
    }
    
    std::vector<DynamicGenomeData>& getGenerationBeforeLastGenomeData(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last data: currentGen underflow");
        return getGenomeData(currentGen - 2);
    }
    
    std::multimap<FitnessResultType, size_t>& getCurrentFitnessResults(uint32_t currentGen) {
        return getFitnessResults(currentGen);
    }
    
    std::multimap<FitnessResultType, size_t>& getLastFitnessResults(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation fitness: currentGen underflow");
        return getFitnessResults(currentGen - 1);
    }
    
    std::multimap<FitnessResultType, size_t>& getGenerationBeforeLastFitnessResults(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last fitness: currentGen underflow");
        return getFitnessResults(currentGen - 2);
    }
    
    // Synchronized addition to maintain vector consistency and registry synchronization
    uint32_t push_back(uint32_t generation, Genome genome, DynamicGenomeData data) {
        size_t targetGeneration = generation % 3;
        
        // Grow registry to match new size
        uint32_t globalIndex = _registry.incrementMaxIndex();
        
        // Add to target generation
        _genomes[targetGeneration].emplace_back(std::move(genome));
        _genomeData[targetGeneration].emplace_back(std::move(data));
        
        // Ensure all other generations have matching sizes by adding copies
        size_t targetSize = _genomes[targetGeneration].size();
        for (size_t i = 0; i < 3; ++i) {
            if (i != targetGeneration) {
                while (_genomes[i].size() < targetSize) {
                    _genomes[i].emplace_back(_genomes[targetGeneration].back());
                    _genomeData[i].emplace_back(_genomeData[targetGeneration].back());
                }
            }
        }
        
        return globalIndex;
    }
    
    // Strict replacement at specific index with validation
    void replace(uint32_t generation, size_t targetIndex, Genome genome, DynamicGenomeData data) {
        size_t targetGeneration = generation % 3;
        
        // Strict bounds checking
        assert(targetIndex < _genomes[targetGeneration].size() && 
               "replace() targetIndex out of bounds for target generation");
        assert(targetIndex < _genomeData[targetGeneration].size() && 
               "replace() targetIndex out of bounds for genome data");
        assert(_genomes[targetGeneration].size() == _genomeData[targetGeneration].size() && 
               "Genome and data vectors must have equal size");
        
        // Verify the target index is ready for replacement (Active means it was reserved by getFreeIndex)
        assert((_registry.getState(static_cast<uint32_t>(targetIndex)) == GenomeState::ReadyForReplacement ||
                _registry.getState(static_cast<uint32_t>(targetIndex)) == GenomeState::Active) &&
               "replace() targetIndex must be in ReadyForReplacement or Active (reserved) state");
        
        // Replace at specific index
        _genomes[targetGeneration][targetIndex] = std::move(genome);
        _genomeData[targetGeneration][targetIndex] = std::move(data);
        
        // Propagate to other generations to maintain consistency
        for (size_t i = 0; i < 3; ++i) {
            if (i != targetGeneration) {
                assert(targetIndex < _genomes[i].size() && 
                       "replace() targetIndex out of bounds for generation sync");
                _genomes[i][targetIndex] = _genomes[targetGeneration][targetIndex];
                _genomeData[i][targetIndex] = _genomeData[targetGeneration][targetIndex];
            }
        }
        
        // Ensure registry state is Active after successful replacement
        if (_registry.getState(static_cast<uint32_t>(targetIndex)) != GenomeState::Active) {
            _registry._states[targetIndex] = GenomeState::Active;
        }
    }

private:
    // Triple-buffer storage arrays
    std::array<std::vector<Genome>, 3> _genomes;
    std::array<std::vector<DynamicGenomeData>, 3> _genomeData;
    std::array<std::multimap<FitnessResultType, size_t>, 3> _fitnessResults;
    
    // Reference to GlobalIndexRegistry for synchronization
    GlobalIndexRegistry& _registry;
};