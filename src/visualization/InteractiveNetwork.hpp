// interactive_network.hpp
#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <string>

namespace neat {
namespace core {
class Genome;
}
namespace visualization {
class InteractiveNetwork {
public:
    struct Config {
        int updateInterval = 100;  // milliseconds
        bool showNodeValues = true;
        bool animateSignals = true;
        bool showStatistics = true;
    };

    InteractiveNetwork(const core::Genome& genome, const Config& config = Config{});
    
    // Start interactive visualization
    void start();
    
    // Update network state
    void updateInputs(const std::vector<double>& inputs);
    void step();
    void reset();
    
    // Event handlers
    void onNodeClick(int32_t nodeId);
    void onConnectionClick(int32_t fromId, int32_t toId);
    
    // Component generation
    std::string generateReactComponent();

private:
    const core::Genome& genome;
    Config config;
    std::unordered_map<int32_t, double> nodeValues;
    std::vector<std::vector<double>> signalHistory;
};

}
}
