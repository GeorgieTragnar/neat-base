#include "SVGGenerator.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Visualization {

// Global configuration
static VisualizationConfig g_config;
static std::string g_outputPath;
static bool g_initialized = false;

std::string getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    
    std::stringstream ss;
    ss << std::put_time(&tm, "%y_%m_%d_%H_%M_%S");
    return ss.str();
}

std::string createOutputDirectory() {
    std::string timestamp = getCurrentTimestamp();
    std::filesystem::path outputDir = std::filesystem::current_path() / g_config.outputDirectory / timestamp;
    
    try {
        std::filesystem::create_directories(outputDir);
    } catch (const std::filesystem::filesystem_error& e) {
        throw std::runtime_error("Failed to create output directory: " + std::string(e.what()));
    }
    
    return outputDir.string();
}

std::string generateFilename(size_t generation, size_t speciesIndex) {
    std::stringstream ss;
    ss << "generation_" << generation << "_species_" << speciesIndex << ".html";
    return ss.str();
}

void initialize(const VisualizationConfig& config) {
    g_config = config;
    g_outputPath = createOutputDirectory();
    g_initialized = true;
}

// Layout structures
struct NodePosition {
    size_t nodeIndex;
    double x, y;
    NodeType type;
    NodeGeneAttributes attributes;
};

struct ConnectionLayout {
    size_t sourceIndex, targetIndex;
    NodePosition source, target;
    Phenotype::Connection connection;
};

class LayoutEngine {
public:
    struct LayoutResult {
        std::vector<NodePosition> nodes;
        std::vector<ConnectionLayout> connections;
        double totalWidth, totalHeight;
    };
    
    LayoutResult computeLayout(const Phenotype& phenotype) {
        LayoutResult result;
        
        if (phenotype._nodeGeneAttributes.empty()) {
            result.totalWidth = g_config.canvasWidth;
            result.totalHeight = g_config.canvasHeight;
            return result;
        }
        
        // Assign nodes to layers
        auto layers = assignLayers(phenotype);
        
        // Calculate positions
        calculatePositions(layers, phenotype, result);
        
        // Create connection layouts
        createConnectionLayouts(phenotype, result);
        
        return result;
    }
    
private:
    std::vector<std::vector<size_t>> assignLayers(const Phenotype& phenotype) {
        const auto& inputIndices = phenotype._inputIndices;
        const auto& outputIndices = phenotype._outputIndices;
        
        std::vector<std::vector<size_t>> layers;
        
        // Layer 0: Inputs and bias (non-output nodes)
        std::vector<size_t> inputLayer;
        for (size_t i = 0; i < phenotype._nodeGeneAttributes.size(); ++i) {
            if (std::find(outputIndices.begin(), outputIndices.end(), i) == outputIndices.end()) {
                inputLayer.push_back(i);
            }
        }
        layers.push_back(inputLayer);
        
        // Final layer: Outputs
        layers.push_back(outputIndices);
        
        return layers;
    }
    
    void calculatePositions(const std::vector<std::vector<size_t>>& layers, 
                           const Phenotype& phenotype,
                           LayoutResult& result) {
        
        // Calculate canvas dimensions based on layers
        size_t maxLayerSize = 0;
        for (const auto& layer : layers) {
            maxLayerSize = std::max(maxLayerSize, layer.size());
        }
        
        double totalWidth = layers.size() * g_config.layerSpacing;
        double totalHeight = maxLayerSize * g_config.nodeSpacing + 2 * g_config.nodeRadius;
        
        result.totalWidth = std::max(totalWidth, g_config.canvasWidth);
        result.totalHeight = std::max(totalHeight, g_config.canvasHeight);
        
        // Position nodes
        for (size_t layerIdx = 0; layerIdx < layers.size(); ++layerIdx) {
            const auto& layer = layers[layerIdx];
            double x = g_config.nodeRadius + layerIdx * g_config.layerSpacing;
            
            // Center the layer vertically
            double layerHeight = layer.size() * g_config.nodeSpacing;
            double startY = (result.totalHeight - layerHeight) / 2.0 + g_config.nodeRadius;
            
            for (size_t nodeIdx = 0; nodeIdx < layer.size(); ++nodeIdx) {
                size_t globalNodeIdx = layer[nodeIdx];
                double y = startY + nodeIdx * g_config.nodeSpacing;
                
                NodePosition pos;
                pos.nodeIndex = globalNodeIdx;
                pos.x = x;
                pos.y = y;
                pos.attributes = phenotype._nodeGeneAttributes[globalNodeIdx];
                
                // Determine node type based on indices
                const auto& inputIndices = phenotype._inputIndices;
                const auto& outputIndices = phenotype._outputIndices;
                
                if (std::find(inputIndices.begin(), inputIndices.end(), globalNodeIdx) != inputIndices.end()) {
                    pos.type = NodeType::INPUT;
                } else if (std::find(outputIndices.begin(), outputIndices.end(), globalNodeIdx) != outputIndices.end()) {
                    pos.type = NodeType::OUTPUT;
                } else {
                    // For now, assume all non-input, non-output nodes are bias
                    // This is a simplification - in full NEAT, we'd need to track node types
                    pos.type = NodeType::BIAS;
                }
                
                result.nodes.push_back(pos);
            }
        }
    }
    
    void createConnectionLayouts(const Phenotype& phenotype, LayoutResult& result) {
        for (const auto& conn : phenotype._orderedConnections) {
            if (conn._sourceNodeIndex < result.nodes.size() && 
                conn._targetNodeIndex < result.nodes.size()) {
                
                ConnectionLayout layout;
                layout.sourceIndex = conn._sourceNodeIndex;
                layout.targetIndex = conn._targetNodeIndex;
                layout.source = result.nodes[conn._sourceNodeIndex];
                layout.target = result.nodes[conn._targetNodeIndex];
                layout.connection = conn;
                
                result.connections.push_back(layout);
            }
        }
    }
};

class SVGBuilder {
public:
    void addNode(const NodePosition& node) {
        std::string color = getNodeColor(node.type);
        std::string shape = (node.type == NodeType::BIAS) ? "rect" : "circle";
        
        if (shape == "circle") {
            svgContent << "<circle cx=\"" << node.x << "\" cy=\"" << node.y 
                      << "\" r=\"" << g_config.nodeRadius << "\" fill=\"" << color 
                      << "\" stroke=\"#000\" stroke-width=\"1\"";
        } else {
            double size = g_config.nodeRadius * 1.5;
            svgContent << "<rect x=\"" << (node.x - size/2) << "\" y=\"" << (node.y - size/2)
                      << "\" width=\"" << size << "\" height=\"" << size 
                      << "\" fill=\"" << color << "\" stroke=\"#000\" stroke-width=\"1\"";
        }
        
        // Add tooltip
        svgContent << " title=\"";
        if (node.type == NodeType::BIAS) {
            svgContent << "Node Type: BIAS";
        } else {
            svgContent << "ActivationType: " << activationTypeToString(node.attributes.activationType);
        }
        svgContent << "\"/>\n";
    }
    
    void addConnection(const ConnectionLayout& conn) {
        double weight = conn.connection._connectionGeneAttribute.weight;
        bool enabled = conn.connection._connectionGeneAttribute.enabled;
        
        double lineWidth = calculateLineWidth(weight);
        std::string color = enabled ? g_config.enabledLineColor : g_config.disabledLineColor;
        std::string dashArray = enabled ? "" : " stroke-dasharray=\"5,5\"";
        
        // Calculate arrow position
        double dx = conn.target.x - conn.source.x;
        double dy = conn.target.y - conn.source.y;
        double length = std::sqrt(dx*dx + dy*dy);
        
        if (length > 0) {
            dx /= length;
            dy /= length;
            
            // Adjust start and end points to node boundaries
            double startX = conn.source.x + dx * g_config.nodeRadius;
            double startY = conn.source.y + dy * g_config.nodeRadius;
            double endX = conn.target.x - dx * g_config.nodeRadius;
            double endY = conn.target.y - dy * g_config.nodeRadius;
            
            svgContent << "<line x1=\"" << startX << "\" y1=\"" << startY
                      << "\" x2=\"" << endX << "\" y2=\"" << endY
                      << "\" stroke=\"" << color << "\" stroke-width=\"" << lineWidth << "\""
                      << dashArray;
            
            // Add tooltip with weight
            svgContent << " title=\"Weight: " << std::fixed << std::setprecision(10) << weight << "\"/>\n";
            
            // Add arrowhead
            addArrowhead(endX, endY, dx, dy, color);
        }
    }
    
    std::string generateSVG(double width, double height) {
        std::stringstream svg;
        svg << "<svg width=\"" << width << "\" height=\"" << height 
            << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        svg << svgContent.str();
        svg << "</svg>\n";
        return svg.str();
    }
    
private:
    std::stringstream svgContent;
    
    std::string getNodeColor(NodeType type) {
        switch (type) {
            case NodeType::INPUT: return g_config.inputNodeColor;
            case NodeType::OUTPUT: return g_config.outputNodeColor;
            case NodeType::BIAS: return g_config.biasNodeColor;
            case NodeType::HIDDEN: return g_config.hiddenNodeColor;
            default: return g_config.hiddenNodeColor;
        }
    }
    
    std::string activationTypeToString(ActivationType type) {
        switch (type) {
            case ActivationType::NONE: return "NONE";
            case ActivationType::SIGMOID: return "SIGMOID";
            case ActivationType::TANH: return "TANH";
            case ActivationType::RELU: return "RELU";
            case ActivationType::LEAKY_RELU: return "LEAKY_RELU";
            case ActivationType::STEP: return "STEP";
            default: return "UNKNOWN";
        }
    }
    
    double calculateLineWidth(float weight) {
        double magnitude = std::abs(weight);
        double normalized = std::min(1.0, magnitude / 2.0); // Normalize to 0-1 range
        return g_config.minLineWidth + normalized * (g_config.maxLineWidth - g_config.minLineWidth);
    }
    
    void addArrowhead(double x, double y, double dx, double dy, const std::string& color) {
        double arrowSize = 8.0;
        double angle = std::atan2(dy, dx);
        
        double x1 = x - arrowSize * std::cos(angle - 0.5);
        double y1 = y - arrowSize * std::sin(angle - 0.5);
        double x2 = x - arrowSize * std::cos(angle + 0.5);
        double y2 = y - arrowSize * std::sin(angle + 0.5);
        
        svgContent << "<polygon points=\"" << x << "," << y << " " 
                   << x1 << "," << y1 << " " << x2 << "," << y2 
                   << "\" fill=\"" << color << "\"/>\n";
    }
};

void generateVisualization(const Phenotype& phenotype, 
                          size_t generation, 
                          size_t speciesIndex) {
    if (!g_initialized) {
        throw std::runtime_error("Visualization system not initialized. Call initialize() first.");
    }
    
    // Compute layout
    LayoutEngine layoutEngine;
    auto layout = layoutEngine.computeLayout(phenotype);
    
    // Generate SVG
    SVGBuilder svgBuilder;
    
    // Add all nodes
    for (const auto& node : layout.nodes) {
        svgBuilder.addNode(node);
    }
    
    // Add all connections
    for (const auto& conn : layout.connections) {
        svgBuilder.addConnection(conn);
    }
    
    std::string svgContent = svgBuilder.generateSVG(layout.totalWidth, layout.totalHeight);
    
    // Create HTML wrapper
    std::stringstream html;
    html << "<!DOCTYPE html>\n<html>\n<head>\n";
    html << "<title>Generation " << generation << " Species " << speciesIndex << "</title>\n";
    html << "<style>\n";
    html << "body { font-family: Arial, sans-serif; margin: 20px; }\n";
    html << "h1 { color: #333; }\n";
    html << "svg { border: 1px solid #ccc; }\n";
    html << "</style>\n";
    html << "</head>\n<body>\n";
    html << "<h1>Phenotype Visualization - Generation " << generation << ", Species " << speciesIndex << "</h1>\n";
    html << svgContent;
    html << "</body>\n</html>\n";
    
    // Write to file
    std::string filename = generateFilename(generation, speciesIndex);
    std::filesystem::path filepath = std::filesystem::path(g_outputPath) / filename;
    
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to create output file: " + filepath.string());
    }
    
    file << html.str();
    file.close();
}

}