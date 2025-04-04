// InnovationTracker.cpp
#include "InnovationTracker.hpp"

namespace neat {
namespace core {
    int32_t InnovationTracker::nextInnovation = 0;
    std::unordered_map<
        std::pair<int32_t, int32_t>, 
        int32_t, 
        InnovationTracker::PairHash
    > InnovationTracker::currentGenerationInnovations;
}
}
