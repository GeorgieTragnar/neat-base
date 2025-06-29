#pragma once

#include <array>
#include <vector>
#include <map>
#include <cstdint>
#include <cassert>

// Forward declarations
class Genome;

namespace Population {

// Forward declaration
struct DynamicGenomeData;

template<typename FitnessResultType>
class PopulationContainer {
public:
    PopulationContainer() = default;
    ~PopulationContainer() = default;

    // Core generation-based access using modulo arithmetic
    std::vector<Genome>& getGenomes(uint32_t generation) {
        return _genomes[generation % 3];
    }
    
    const std::vector<Genome>& getGenomes(uint32_t generation) const {
        return _genomes[generation % 3];
    }
    
    std::vector<DynamicGenomeData>& getGenomeData(uint32_t generation) {
        return _genomeData[generation % 3];
    }
    
    const std::vector<DynamicGenomeData>& getGenomeData(uint32_t generation) const {
        return _genomeData[generation % 3];
    }
    
    std::multimap<FitnessResultType, size_t>& getFitnessResults(uint32_t generation) {
        return _fitnessResults[generation % 3];
    }
    
    const std::multimap<FitnessResultType, size_t>& getFitnessResults(uint32_t generation) const {
        return _fitnessResults[generation % 3];
    }
    
    // Convenience methods for common access patterns with underflow protection
    std::vector<Genome>& getCurrentGenomes(uint32_t currentGen) {
        return getGenomes(currentGen);
    }
    
    const std::vector<Genome>& getCurrentGenomes(uint32_t currentGen) const {
        return getGenomes(currentGen);
    }
    
    std::vector<Genome>& getLastGenomes(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation: currentGen underflow");
        return getGenomes(currentGen - 1);
    }
    
    const std::vector<Genome>& getLastGenomes(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation: currentGen underflow");
        return getGenomes(currentGen - 1);
    }
    
    std::vector<Genome>& getGenerationBeforeLastGenomes(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last: currentGen underflow");
        return getGenomes(currentGen - 2);
    }
    
    const std::vector<Genome>& getGenerationBeforeLastGenomes(uint32_t currentGen) const {
        assert(currentGen > 1 && "Cannot access generation before last: currentGen underflow");
        return getGenomes(currentGen - 2);
    }
    
    // Genome data convenience methods with underflow protection
    std::vector<DynamicGenomeData>& getCurrentGenomeData(uint32_t currentGen) {
        return getGenomeData(currentGen);
    }
    
    const std::vector<DynamicGenomeData>& getCurrentGenomeData(uint32_t currentGen) const {
        return getGenomeData(currentGen);
    }
    
    std::vector<DynamicGenomeData>& getLastGenomeData(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation data: currentGen underflow");
        return getGenomeData(currentGen - 1);
    }
    
    const std::vector<DynamicGenomeData>& getLastGenomeData(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation data: currentGen underflow");
        return getGenomeData(currentGen - 1);
    }
    
    std::vector<DynamicGenomeData>& getGenerationBeforeLastGenomeData(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last data: currentGen underflow");
        return getGenomeData(currentGen - 2);
    }
    
    const std::vector<DynamicGenomeData>& getGenerationBeforeLastGenomeData(uint32_t currentGen) const {
        assert(currentGen > 1 && "Cannot access generation before last data: currentGen underflow");
        return getGenomeData(currentGen - 2);
    }
    
    // Fitness results convenience methods with underflow protection
    std::multimap<FitnessResultType, size_t>& getCurrentFitnessResults(uint32_t currentGen) {
        return getFitnessResults(currentGen);
    }
    
    const std::multimap<FitnessResultType, size_t>& getCurrentFitnessResults(uint32_t currentGen) const {
        return getFitnessResults(currentGen);
    }
    
    std::multimap<FitnessResultType, size_t>& getLastFitnessResults(uint32_t currentGen) {
        assert(currentGen > 0 && "Cannot access last generation fitness: currentGen underflow");
        return getFitnessResults(currentGen - 1);
    }
    
    const std::multimap<FitnessResultType, size_t>& getLastFitnessResults(uint32_t currentGen) const {
        assert(currentGen > 0 && "Cannot access last generation fitness: currentGen underflow");
        return getFitnessResults(currentGen - 1);
    }
    
    std::multimap<FitnessResultType, size_t>& getGenerationBeforeLastFitnessResults(uint32_t currentGen) {
        assert(currentGen > 1 && "Cannot access generation before last fitness: currentGen underflow");
        return getFitnessResults(currentGen - 2);
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

private:
    // Triple-buffer storage arrays
    std::array<std::vector<Genome>, 3> _genomes;
    std::array<std::vector<DynamicGenomeData>, 3> _genomeData;
    std::array<std::multimap<FitnessResultType, size_t>, 3> _fitnessResults;
};

} // namespace Population