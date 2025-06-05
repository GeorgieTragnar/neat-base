#pragma once
#include "FitnessResult.hpp"
#include "SpeciationControlUnit.hpp"

namespace Analysis {

template<typename FitnessResultType>
class FitnessStrategy {
public:
    virtual ~FitnessStrategy() = default;
    
    virtual FitnessResultType evaluate(
        const Phenotype& candidate,
        const SpeciationControlUnit& speciation) const = 0;
};

}