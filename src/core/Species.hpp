// species.hpp
#pragma once
#include <cstdint>
#include <vector>
#include "Gene.hpp"

namespace neat {
namespace core {

class Species {
public:
    using GenomeIndex = int32_t;
    using Fitness = double;
    
    void addGenome(GenomeIndex idx);
    void removeGenome(GenomeIndex idx);
    void updateRepresentative(const class Genome& genome);
    void updateStats();
    
    Fitness getBestFitness() const noexcept { return bestFitness; }
    int32_t getStaleness() const noexcept { return staleness; }
    const std::vector<GenomeIndex>& getMembers() const noexcept { return members; }

private:
    std::vector<GenomeIndex> members;
    Fitness bestFitness = 0.0;
    int32_t staleness = 0;
    Gene representative;
    double totalAdjustedFitness = 0.0;
};

}
}