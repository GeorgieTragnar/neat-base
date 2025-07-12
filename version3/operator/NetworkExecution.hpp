#pragma once

#include <vector>
#include "version3/data/Phenotype.hpp"

namespace Operator {

// Parameters for network execution
class NetworkExecutionParams {
public:
    NetworkExecutionParams() = delete;
    NetworkExecutionParams(bool debugOutput);

protected:
    friend std::vector<double> networkExecution(
        const Phenotype& phenotype,
        const std::vector<double>& inputs,
        const NetworkExecutionParams& params
    );
    
    const bool _debugOutput;
};

// Execute network and return all output values
std::vector<double> networkExecution(
    const Phenotype& phenotype,
    const std::vector<double>& inputs,
    const NetworkExecutionParams& params
);

} // namespace Operator