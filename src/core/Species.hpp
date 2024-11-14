// species.hpp
#pragma once
#include <cstdint>
#include <vector>
#include <algorithm>
#include "Gene.hpp"

namespace neat {
namespace core {

class Species {
public:
    using GenomeIndex = int32_t;
    using Fitness = double;

    Species() = default;
    
    void addGenome(GenomeIndex idx) { members.push_back(idx); }
    void removeGenome(GenomeIndex idx) {
        members.erase(std::remove(members.begin(), members.end(), idx), members.end());
    }
    void setRepresentative(GenomeIndex idx) { representative = idx; }

    void clearMembers() {
        members.clear();
    }

    void setBestFitness(Fitness f) {
        if (f > bestFitness) {
            bestFitness = f;
            staleness = 0;
        }
    }

    void updateFitness(Fitness newFitness) {
        if (newFitness > bestFitness) {
            bestFitness = newFitness;
            staleness = 0;
        } else {
            staleness++;
        }
    }

    void resetStaleness() { staleness = 0; }
    void incrementStaleness() { ++staleness; }
    
    const std::vector<GenomeIndex>& getMembers() const noexcept { return members; }
    std::vector<GenomeIndex>& getMembers() noexcept { return members; }
    GenomeIndex getRepresentative() const { return representative; }
    Fitness getBestFitness() const noexcept { return bestFitness; }
    int32_t getStaleness() const noexcept { return staleness; }

private:
    std::vector<GenomeIndex> members;
    Fitness bestFitness = 0.0;
    int32_t staleness = 0;
    GenomeIndex representative;
};

}
}