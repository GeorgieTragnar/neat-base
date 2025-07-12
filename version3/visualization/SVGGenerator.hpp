#pragma once

#include <string>
#include <memory>
#include "../data/Phenotype.hpp"

namespace Visualization {

struct VisualizationConfig {
    std::string outputDirectory = "visualizations";
    double nodeRadius = 15.0;
    double layerSpacing = 120.0;
    double nodeSpacing = 40.0;
    double canvasWidth = 800.0;
    double canvasHeight = 600.0;
    
    // Node colors
    std::string inputNodeColor = "#4A90E2";      // Blue
    std::string outputNodeColor = "#7ED321";     // Green  
    std::string biasNodeColor = "#F5A623";       // Orange
    std::string hiddenNodeColor = "#9B9B9B";     // Grey
    
    // Connection styling
    double minLineWidth = 1.0;
    double maxLineWidth = 5.0;
    std::string enabledLineColor = "#333333";
    std::string disabledLineColor = "#CCCCCC";
};

// Main API functions
void initialize(const VisualizationConfig& config = {});
void generateVisualization(const Phenotype& phenotype, 
                          size_t generation, 
                          size_t speciesIndex);

// Utility functions
std::string createOutputDirectory();
std::string generateFilename(size_t generation, size_t speciesIndex);
std::string getCurrentTimestamp();

}