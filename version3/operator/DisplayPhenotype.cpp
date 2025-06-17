#include "DisplayPhenotype.hpp"
#include <iomanip>

namespace Operator {

void displayPhenotype(const Genome& genome, std::stringstream& output) {
    const auto& phenotype = genome.get_phenotype();
    
    output << "=== PHENOTYPE DISPLAY ===\n";
    
    // Display Phenotype Nodes
    output << "\nPhenotype Nodes (" << phenotype._nodeGeneAttributes.size() << "):\n";
    output << std::setw(8) << "Index" << std::setw(12) << "Activation" << "\n";
    output << std::string(20, '-') << "\n";
    
    for (size_t i = 0; i < phenotype._nodeGeneAttributes.size(); ++i) {
        output << std::setw(8) << i;
        
        // Convert ActivationType to string
        switch (phenotype._nodeGeneAttributes[i].activationType) {
            case ActivationType::NONE: output << std::setw(12) << "NONE"; break;
            case ActivationType::SIGMOID: output << std::setw(12) << "SIGMOID"; break;
            case ActivationType::TANH: output << std::setw(12) << "TANH"; break;
            case ActivationType::RELU: output << std::setw(12) << "RELU"; break;
            case ActivationType::LEAKY_RELU: output << std::setw(12) << "LEAKY_RELU"; break;
            case ActivationType::STEP: output << std::setw(12) << "STEP"; break;
            default: output << std::setw(12) << "UNKNOWN"; break;
        }
        
        output << "\n";
    }
    
    // Display Input Indices
    output << "\nInput Indices (" << phenotype._inputIndices.size() << "): ";
    for (size_t i = 0; i < phenotype._inputIndices.size(); ++i) {
        if (i > 0) output << ", ";
        output << phenotype._inputIndices[i];
    }
    output << "\n";
    
    // Display Output Indices
    output << "\nOutput Indices (" << phenotype._outputIndices.size() << "): ";
    for (size_t i = 0; i < phenotype._outputIndices.size(); ++i) {
        if (i > 0) output << ", ";
        output << phenotype._outputIndices[i];
    }
    output << "\n";
    
    // Display Ordered Connections
    output << "\nOrdered Connections (" << phenotype._orderedConnections.size() << "):\n";
    output << std::setw(8) << "Source" << std::setw(8) << "Target" 
           << std::setw(10) << "Weight" << std::setw(10) << "Enabled" << "\n";
    output << std::string(36, '-') << "\n";
    
    for (const auto& conn : phenotype._orderedConnections) {
        output << std::setw(8) << conn._sourceNodeIndex
               << std::setw(8) << conn._targetNodeIndex
               << std::setw(10) << std::fixed << std::setprecision(3) << conn._connectionGeneAttribute.weight
               << std::setw(10) << (conn._connectionGeneAttribute.enabled ? "YES" : "NO")
               << "\n";
    }
    
    output << "\n";
}

}