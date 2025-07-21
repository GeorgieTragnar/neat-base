#pragma once
#include "../data/Genome.hpp"
#include <memory>

namespace Operator {

class CycleDetectionParams;

bool hasCycles(const Genome& genome, const CycleDetectionParams& params);

class CycleDetectionParams {
public:
    // Constructor with default parameters
    CycleDetectionParams() = delete;
    
    // Constructor with explicit configuration options
    CycleDetectionParams(bool includeAllNodeTypes, bool failFast);

protected:
    friend bool hasCycles(const Genome& genome, const CycleDetectionParams& params);
    
    // Whether to include all node types in cycle detection or only HIDDEN nodes
    // Default: true (include INPUT, OUTPUT, BIAS, HIDDEN - full network analysis)
    const bool _includeAllNodeTypes = true;
    
    // Whether to return immediately on first cycle found (fail-fast optimization)
    // Default: true (return as soon as any cycle is detected)
    const bool _failFast = true;
};

} // namespace Operator