// network_renderer.hpp
#pragma once
#include <string>
#include <sstream>
#include <unordered_map>
#include <memory>

namespace neat {
namespace core {
class Genome;
class Gene;
enum class ENodeType;
}
namespace visualization {

class NetworkRenderer {
public:
    struct Config {
        // SVG styling
        struct Style {
            std::string nodeColor = "#2196F3";
            std::string edgeColor = "#757575";
            std::string positiveWeightColor = "#4CAF50";
            std::string negativeWeightColor = "#F44336";
            std::string nodeTextColor = "#FFFFFF";
            std::string backgroundColor = "#FFFFFF";
            
            double nodeRadius = 20.0;
            double edgeWidth = 2.0;
            double arrowSize = 10.0;
            
            bool showWeights = true;
            bool showNodeIds = true;
            bool showBiasNode = true;
            bool colorizeWeights = true;
        };
        
        int width = 800;
        int height = 600;
        Style style;
    };

    explicit NetworkRenderer(const Config& config);
    
    // Generate SVG representation of the network
    std::string renderToSVG(const core::Genome& genome);
    
    // Generate DOT format for use with Graphviz
    std::string renderToDOT(const core::Genome& genome);
    
    // Save to file
    void saveToFile(const core::Genome& genome, const std::string& filename);

private:
    struct NodeLayout {
        double x, y;
        core::ENodeType type;
        int32_t id;
    };

    void calculateLayout(const core::Genome& genome);
    std::string generateNodeSVG(const NodeLayout& node);
    std::string generateEdgeSVG(const NodeLayout& from, const NodeLayout& to, const core::Gene& gene);
    std::string getColorForWeight(double weight);

    Config config;
    std::unordered_map<int32_t, NodeLayout> nodePositions;
};

}
}
