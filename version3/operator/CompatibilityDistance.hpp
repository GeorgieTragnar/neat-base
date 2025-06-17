#pragma once

#include <memory>

#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"

namespace Operator {

class CompatibilityDistanceParams;

uint32_t compatibilityDistance(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const CompatibilityDistanceParams& params);

class CompatibilityDistanceParams {
public:
    CompatibilityDistanceParams() = delete;
    CompatibilityDistanceParams(float c1, float c2, float c3, float threshold);

protected:
    friend uint32_t compatibilityDistance(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const CompatibilityDistanceParams& params);
    
    const float _c1;          // Coefficient for excess genes
    const float _c2;          // Coefficient for disjoint genes  
    const float _c3;          // Coefficient for weight differences
    const float _threshold;   // Compatibility threshold
};

}