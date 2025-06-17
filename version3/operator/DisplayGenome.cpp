#include "DisplayGenome.hpp"
#include <iomanip>

namespace Operator {

void displayGenome(const Genome& genome, std::stringstream& output) {
    const auto& nodeGenes = genome.get_nodeGenes();
    const auto& connectionGenes = genome.get_connectionGenes();
    
    output << "=== GENOME DISPLAY ===\n";
    
    // Display Node Genes
    output << "\nNode Genes (" << nodeGenes.size() << "):\n";
    output << std::setw(8) << "HistID" << std::setw(8) << "Type" << std::setw(12) << "Activation" << "\n";
    output << std::string(28, '-') << "\n";
    
    for (const auto& node : nodeGenes) {
        output << std::setw(8) << node.get_historyID();
        
        // Convert NodeType to string
        switch (node.get_type()) {
            case NodeType::INPUT: output << std::setw(8) << "INPUT"; break;
            case NodeType::OUTPUT: output << std::setw(8) << "OUTPUT"; break;
            case NodeType::HIDDEN: output << std::setw(8) << "HIDDEN"; break;
            case NodeType::BIAS: output << std::setw(8) << "BIAS"; break;
            default: output << std::setw(8) << "UNKNOWN"; break;
        }
        
        // Convert ActivationType to string
        switch (node.get_attributes().activationType) {
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
    
    // Display Connection Genes
    output << "\nConnection Genes (" << connectionGenes.size() << "):\n";
    output << std::setw(8) << "HistID" << std::setw(8) << "From" << std::setw(8) << "To" 
           << std::setw(10) << "Weight" << std::setw(10) << "Enabled" << "\n";
    output << std::string(44, '-') << "\n";
    
    for (const auto& conn : connectionGenes) {
        uint32_t sourceHistoryID = nodeGenes[conn.get_sourceNodeIndex()].get_historyID();
        uint32_t targetHistoryID = nodeGenes[conn.get_targetNodeIndex()].get_historyID();
        
        output << std::setw(8) << conn.get_historyID()
               << std::setw(8) << sourceHistoryID
               << std::setw(8) << targetHistoryID
               << std::setw(10) << std::fixed << std::setprecision(3) << conn.get_attributes().weight
               << std::setw(10) << (conn.get_attributes().enabled ? "YES" : "NO")
               << "\n";
    }
    
    // Display Delta Information (const_cast needed since deltas only have non-const getters)
    const auto& connectionDeltas = const_cast<Genome&>(genome).get_connectionGeneDeltas();
    const auto& nodeDeltas = const_cast<Genome&>(genome).get_nodeGeneDeltas();
    
    if (!connectionDeltas.empty() || !nodeDeltas.empty()) {
        output << "\nDelta Information:\n";
        if (!connectionDeltas.empty()) {
            output << "  Connection Deltas (" << connectionDeltas.size() << "): ";
            for (size_t i = 0; i < connectionDeltas.size(); ++i) {
                if (i > 0) output << ", ";
                output << connectionDeltas[i];
            }
            output << "\n";
        }
        if (!nodeDeltas.empty()) {
            output << "  Node Deltas (" << nodeDeltas.size() << "): ";
            for (size_t i = 0; i < nodeDeltas.size(); ++i) {
                if (i > 0) output << ", ";
                output << nodeDeltas[i];
            }
            output << "\n";
        }
    }
    
    output << "\n";
}

}