#pragma once

#include <vector>
#include "version3/data/Phenotype.hpp"

namespace Operator {

class NetworkExecutionParams;

std::vector<double> networkExecution(
    const Phenotype& phenotype,
    const std::vector<double>& inputs,
    const NetworkExecutionParams& params
);

class NetworkExecutionParams {
public:
    NetworkExecutionParams() = delete;
    NetworkExecutionParams(bool debugOutput);

protected:
    friend std::vector<double> networkExecution(const Phenotype& phenotype, const std::vector<double>& inputs, const NetworkExecutionParams& params);
    
    const bool _debugOutput;
};

} // namespace Operator