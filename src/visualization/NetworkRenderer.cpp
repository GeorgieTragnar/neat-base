#include "NetworkRenderer.hpp"
#include "core/Genome.hpp"
#include "core/Gene.hpp"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "logger/Logger.hpp"
static auto logger = LOGGER("visualization::NetworkRenderer");

namespace neat {
namespace visualization {

NetworkRenderer::NetworkRenderer(const Config& config) : config(config) {}

void NetworkRenderer::calculateLayout(const core::Genome& genome) {
    nodePositions.clear();
    
    // First identify actually connected nodes
    std::set<int32_t> connectedNodes;
    for (const auto& gene : genome.getGenes()) {
        if (gene.enabled) {
            connectedNodes.insert(gene.inputNode);
            connectedNodes.insert(gene.outputNode);
        }
    }
    
    // Calculate horizontal positions
    double inputX = config.width * 0.2;    // 20% from left
    double hiddenX = config.width * 0.5;   // Center
    double outputX = config.width * 0.8;   // 80% from left
    
    // Group nodes by type
    std::vector<int32_t> inputNodes, outputNodes, hiddenNodes;
    for (const auto& [id, type] : genome.getNodes()) {
        if (connectedNodes.count(id) || type == core::ENodeType::INPUT || type == core::ENodeType::OUTPUT) {
            switch (type) {
                case core::ENodeType::INPUT:
                    inputNodes.push_back(id);
                    break;
                case core::ENodeType::OUTPUT:
                    outputNodes.push_back(id);
                    break;
                case core::ENodeType::HIDDEN:
                    if (connectedNodes.count(id)) {
                        hiddenNodes.push_back(id);
                    }
                    break;
                default:
                    break;
            }
        }
    }
    
    // Sort nodes by ID for consistent layout
    std::sort(inputNodes.begin(), inputNodes.end());
    std::sort(hiddenNodes.begin(), hiddenNodes.end());
    std::sort(outputNodes.begin(), outputNodes.end());
    
    // Calculate vertical spacing
    double verticalSpacing = config.height / (std::max({inputNodes.size(), hiddenNodes.size(), outputNodes.size()}) + 1);
    
    // Position bias node at top left
    NodeLayout biasLayout;
    biasLayout.id = -1;
    biasLayout.type = core::ENodeType::BIAS;
    biasLayout.x = inputX;
    biasLayout.y = verticalSpacing * 0.5;
    nodePositions[-1] = biasLayout;
    
    // Position input nodes
    for (size_t i = 0; i < inputNodes.size(); ++i) {
        NodeLayout layout;
        layout.id = inputNodes[i];
        layout.type = genome.getNodes().at(inputNodes[i]);
        layout.x = inputX;
        layout.y = (i + 1) * verticalSpacing;
        nodePositions[inputNodes[i]] = layout;
    }
    
    // Position hidden nodes
    for (size_t i = 0; i < hiddenNodes.size(); ++i) {
        NodeLayout layout;
        layout.id = hiddenNodes[i];
        layout.type = genome.getNodes().at(hiddenNodes[i]);
        layout.x = hiddenX;
        layout.y = (i + 1) * verticalSpacing;
        nodePositions[hiddenNodes[i]] = layout;
    }
    
    // Position output nodes
    for (size_t i = 0; i < outputNodes.size(); ++i) {
        NodeLayout layout;
        layout.id = outputNodes[i];
        layout.type = genome.getNodes().at(outputNodes[i]);
        layout.x = outputX;
        layout.y = (i + 1) * verticalSpacing;
        nodePositions[outputNodes[i]] = layout;
    }
    
    LOG_INFO("Layout calculated with:");
    LOG_INFO("  Input nodes: {}", inputNodes.size());
    LOG_INFO("  Hidden nodes (connected): ", hiddenNodes.size());
    LOG_INFO("  Output nodes: {}", outputNodes.size());
}

std::string NetworkRenderer::generateNodeSVG(const NodeLayout& node) {
    std::ostringstream ss;
    double radius = config.style.nodeRadius;
    
    ss << "<circle "
       << "cx=\"" << node.x << "\" "
       << "cy=\"" << node.y << "\" "
       << "r=\"" << radius << "\" "
       << "fill=\"" << config.style.nodeColor << "\" "
       << "/>\n";
       
    if (config.style.showNodeIds) {
        ss << "<text "
           << "x=\"" << node.x << "\" "
           << "y=\"" << node.y << "\" "
           << "text-anchor=\"middle\" "
           << "dominant-baseline=\"middle\" "
           << "fill=\"" << config.style.nodeTextColor << "\" "
           << ">" << node.id << "</text>\n";
    }
    
    return ss.str();
}

std::string NetworkRenderer::generateEdgeSVG(const NodeLayout& from, const NodeLayout& to, const core::Gene& gene) {
    std::ostringstream ss;
    
    // Determine if this is a feedback connection
    bool isFeedback = from.x >= to.x;
    
    // Calculate control points
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    double ctrlX1, ctrlY1, ctrlX2, ctrlY2;
    
    if (isFeedback) {
        // Make feedback connections curve more dramatically
        ctrlX1 = from.x + std::abs(dy) * 0.5;
        ctrlY1 = from.y;
        ctrlX2 = to.x - std::abs(dy) * 0.5;
        ctrlY2 = to.y;
    } else {
        // Forward connections curve gently
        ctrlX1 = from.x + dx * 0.5;
        ctrlY1 = from.y;
        ctrlX2 = to.x - dx * 0.5;
        ctrlY2 = to.y;
    }
    
    // Determine color based on weight
    std::string color = gene.weight > 0 ? config.style.positiveWeightColor : config.style.negativeWeightColor;
    
    // Draw curved connection
    ss << "<path d=\"M " << from.x << " " << from.y << " "
       << "C " << ctrlX1 << " " << ctrlY1 << " "
       << ctrlX2 << " " << ctrlY2 << " "
       << to.x << " " << to.y << "\" "
       << "stroke=\"" << color << "\" "
       << "stroke-width=\"" << config.style.edgeWidth << "\" "
       << "fill=\"none\" "
       << (gene.enabled ? "" : "stroke-dasharray=\"5,5\" ")
       << "/>\n";
    
    // Add weight label
    if (config.style.showWeights) {
        double midX = (from.x + to.x) / 2;
        double midY = (from.y + to.y) / 2;
        if (isFeedback) {
            // Adjust label position for feedback connections
            midX += std::abs(dy) * 0.25;
        }
        ss << "<text "
           << "x=\"" << midX << "\" "
           << "y=\"" << midY << "\" "
           << "text-anchor=\"middle\" "
           << "dominant-baseline=\"middle\" "
           << "fill=\"" << color << "\" "
           << ">" << std::fixed << std::setprecision(2) << gene.weight << "</text>\n";
    }
    
    return ss.str();
}

std::string NetworkRenderer::getColorForWeight(double weight) {
    if (weight > 0) {
        return config.style.positiveWeightColor;
    } else {
        return config.style.negativeWeightColor;
    }
}

std::string NetworkRenderer::renderToSVG(const core::Genome& genome) {
    calculateLayout(genome);
    
    std::ostringstream ss;
    ss << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
       << "width=\"" << config.width << "\" height=\"" << config.height << "\" "
       << "viewBox=\"0 0 " << config.width << " " << config.height << "\">\n";
       
    // Background
    ss << "<rect width=\"100%\" height=\"100%\" fill=\"" 
       << config.style.backgroundColor << "\"/>\n";
    
    // Draw connections
    for (const auto& gene : genome.getGenes()) {
        auto fromIt = nodePositions.find(gene.inputNode);
        auto toIt = nodePositions.find(gene.outputNode);
        if (fromIt != nodePositions.end() && toIt != nodePositions.end()) {
            ss << generateEdgeSVG(fromIt->second, toIt->second, gene);
        }
    }
    
    // Draw nodes
    for (const auto& [id, layout] : nodePositions) {
        if (layout.type != core::ENodeType::BIAS || config.style.showBiasNode) {
            ss << generateNodeSVG(layout);
        }
    }
    
    ss << "</svg>";
    return ss.str();
}

std::string NetworkRenderer::renderToDOT(const core::Genome& genome) {
    std::ostringstream ss;
    ss << "digraph neural_network {\n";
    
    // Node definitions
    for (const auto& [id, type] : genome.getNodes()) {
        ss << "  " << id << " [label=\"" << id << "\"";
        switch (type) {
            case core::ENodeType::INPUT:
                ss << ", shape=circle, style=filled, fillcolor=lightblue";
                break;
            case core::ENodeType::HIDDEN:
                ss << ", shape=circle, style=filled, fillcolor=gray";
                break;
            case core::ENodeType::OUTPUT:
                ss << ", shape=circle, style=filled, fillcolor=lightgreen";
                break;
            case core::ENodeType::BIAS:
                ss << ", shape=diamond, style=filled, fillcolor=yellow";
                break;
        }
        ss << "];\n";
    }
    
    // Edge definitions
    for (const auto& gene : genome.getGenes()) {
        if (gene.enabled) {
            ss << "  " << gene.inputNode << " -> " << gene.outputNode
               << " [label=\"" << std::fixed << std::setprecision(2) << gene.weight << "\"";
            if (config.style.colorizeWeights) {
                ss << ", color=\"" << getColorForWeight(gene.weight) << "\"";
            }
            ss << "];\n";
        }
    }
    
    ss << "}\n";
    return ss.str();
}

void NetworkRenderer::saveToFile(const core::Genome& genome, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << renderToSVG(genome);
        file.close();
    }
}

}
}
