#pragma once

#include "FitnessResult.hpp"

namespace Analysis {

class GenomeAnalytics {
public:
    GenomeAnalytics(const FitnessResult& fitness) : _fitness(fitness) {}
    
    const FitnessResult& getFitness() const {
        return _fitness;
    }

private:
    FitnessResult _fitness;
};

}