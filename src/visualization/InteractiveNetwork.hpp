// interactive_network.hpp
#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <string>
#include "network/Network.hpp"

namespace neat {
namespace core {
class Genome;
}
namespace visualization {
class InteractiveNetwork {
public:
    struct Config {
        int updateInterval;
        bool showNodeValues;
        bool animateSignals;
        bool showStatistics;
        
        Config() 
            : updateInterval(100)
            , showNodeValues(true)
            , animateSignals(true)
            , showStatistics(true) {}
    };

    explicit InteractiveNetwork(const core::Genome& genome, const Config& config = Config{});
    
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
    network::Network::Config netConfig;
    std::unordered_map<int32_t, double> nodeValues;
    std::vector<std::vector<double>> signalHistory;
};

}
}
