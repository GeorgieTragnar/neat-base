#pragma once
#include <functional>
#include <vector>
#include "../FitnessStrategy.hpp"

namespace Analysis {

template<typename FitnessResultType>
struct MultiRunSimpleFitnessParams {
    std::function<FitnessResultType(const Phenotype&)> evaluationFunction;
    int numRuns;
    std::function<FitnessResultType(const std::vector<FitnessResultType>&)> resolver;
};

template<typename FitnessResultType>
class MultiRunSimpleFitnessStrategy : public FitnessStrategy<FitnessResultType> {
public:
    explicit MultiRunSimpleFitnessStrategy(const MultiRunSimpleFitnessParams<FitnessResultType>& params)
        : _params(params) {}
    
    FitnessResultType evaluate(
        const Phenotype& candidate,
        const SpeciationControlUnit& speciation) const override {
        std::vector<FitnessResultType> results;
        results.reserve(_params.numRuns);
        
        for (int i = 0; i < _params.numRuns; ++i) {
            results.push_back(_params.evaluationFunction(candidate));
        }
        
        return _params.resolver(results);
    }

private:
    MultiRunSimpleFitnessParams<FitnessResultType> _params;
};

}