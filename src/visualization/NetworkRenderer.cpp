#include "NetworkRenderer.hpp"
#include "core/Genome.hpp"
#include "core/Gene.hpp"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

namespace neat {
namespace visualization {

NetworkRenderer::NetworkRenderer(const Config& config) : config(config) {}

void NetworkRenderer::calculateLayout(const core::Genome& genome) {
    nodePositions.clear();
    const auto& nodes = genome.getNodes();
    
    // Count nodes of each type
    size_t inputCount = 0, hiddenCount = 0, outputCount = 0;
    for (const auto& [id, type] : nodes) {
        switch (type) {
            case core::ENodeType::INPUT: inputCount++; break;
            case core::ENodeType::HIDDEN: hiddenCount++; break;
            case core::ENodeType::OUTPUT: outputCount++; break;
            default: break;
        }
    }
    
    // Calculate positions
    double verticalSpacing = config.height / std::max((size_t)2, std::max(inputCount, outputCount));
    double horizontalSpacing = config.width / 4.0;
    
    size_t inputIdx = 0, hiddenIdx = 0, outputIdx = 0;
    for (const auto& [id, type] : nodes) {
        NodeLayout layout;
        layout.id = id;
        layout.type = type;
        
        switch (type) {
            case core::ENodeType::INPUT:
                layout.x = horizontalSpacing;
                layout.y = (inputIdx + 1) * verticalSpacing;
                inputIdx++;
                break;
                
            case core::ENodeType::HIDDEN:
                layout.x = 2 * horizontalSpacing;
                layout.y = (hiddenIdx + 1) * (config.height / (hiddenCount + 1));
                hiddenIdx++;
                break;
                
            case core::ENodeType::OUTPUT:
                layout.x = 3 * horizontalSpacing;
                layout.y = (outputIdx + 1) * verticalSpacing;
                outputIdx++;
                break;
                
            case core::ENodeType::BIAS:
                layout.x = horizontalSpacing;
                layout.y = 0;
                break;
        }
        
        nodePositions[id] = layout;
    }
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
    
    // Calculate edge path
    double dx = to.x - from.x;
    double dy = to.y - from.y;
    double angle = std::atan2(dy, dx);
    
    double fromX = from.x + config.style.nodeRadius * std::cos(angle);
    double fromY = from.y + config.style.nodeRadius * std::sin(angle);
    double toX = to.x - config.style.nodeRadius * std::cos(angle);
    double toY = to.y - config.style.nodeRadius * std::sin(angle);
    
    std::string color = config.style.colorizeWeights ? 
        getColorForWeight(gene.weight) : config.style.edgeColor;
    
    // Draw connection line
    ss << "<line "
       << "x1=\"" << fromX << "\" "
       << "y1=\"" << fromY << "\" "
       << "x2=\"" << toX << "\" "
       << "y2=\"" << toY << "\" "
       << "stroke=\"" << color << "\" "
       << "stroke-width=\"" << config.style.edgeWidth << "\" "
       << (gene.enabled ? "" : "stroke-dasharray=\"5,5\" ")
       << "/>\n";
       
    // Add weight label
    if (config.style.showWeights) {
        double midX = (fromX + toX) / 2;
        double midY = (fromY + toY) / 2;
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
