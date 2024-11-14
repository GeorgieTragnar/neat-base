// react_components.hpp
#pragma once

#include <string>

namespace neat {
namespace core {
class Genome;
}
namespace visualization {

class StatsTracker;

struct ReactComponents {
    static std::string generateNetworkViewer(const core::Genome& genome);
    static std::string generateStatsViewer(const StatsTracker& stats);
};

} // namespace visualization
} // namespace neat
