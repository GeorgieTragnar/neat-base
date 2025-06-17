#pragma once
#include <functional>
#include "../FitnessStrategy.hpp"

namespace Analysis {

template<typename FitnessResultType>
struct SingleSimpleFitnessParams {
    std::function<FitnessResultType(const Phenotype&)> evaluationFunction;
};

template<typename FitnessResultType>
class SingleSimpleFitnessStrategy : public FitnessStrategy<FitnessResultType> {
public:
    explicit SingleSimpleFitnessStrategy(const SingleSimpleFitnessParams<FitnessResultType>& params)
        : _params(params) {}
    
    FitnessResultType evaluate(
        const Phenotype& candidate,
        const SpeciationControlUnit& speciation) const override {
        return _params.evaluationFunction(candidate);
    }

private:
    SingleSimpleFitnessParams<FitnessResultType> _params;
};

}