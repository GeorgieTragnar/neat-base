#pragma once
#include "../data/Genome.hpp"
#include <memory>

namespace Operator {

class CycleDetectionParams;

// Cycle detection function - detects cycles in enabled connection graph
bool hasCycles(const Genome& genome, const CycleDetectionParams& params);

class CycleDetectionParams {
public:
    // Constructor with default parameters
    CycleDetectionParams() = default;
    
    // Constructor with explicit configuration options
    explicit CycleDetectionParams(bool includeAllNodeTypes, bool failFast = true);

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