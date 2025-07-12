#include "GenomeEquals.hpp"

namespace Operator {

bool genomeEquals(const Genome& genome1, const Genome& genome2) {
    const auto& nodes1 = genome1.get_nodeGenes();
    const auto& nodes2 = genome2.get_nodeGenes();
    const auto& connections1 = genome1.get_connectionGenes();
    const auto& connections2 = genome2.get_connectionGenes();
    
    // Check node count
    if (nodes1.size() != nodes2.size()) {
        return false;
    }
    
    // Check connection count
    if (connections1.size() != connections2.size()) {
        return false;
    }
    
    // Check all nodes
    for (size_t i = 0; i < nodes1.size(); ++i) {
        if (nodes1[i].get_historyID() != nodes2[i].get_historyID()) {
            return false;
        }
        if (nodes1[i].get_type() != nodes2[i].get_type()) {
            return false;
        }
        if (nodes1[i].get_attributes().activationType != nodes2[i].get_attributes().activationType) {
            return false;
        }
    }
    
    // Check all connections
    for (size_t i = 0; i < connections1.size(); ++i) {
        if (connections1[i].get_historyID() != connections2[i].get_historyID()) {
            return false;
        }
        if (connections1[i].get_attributes().weight != connections2[i].get_attributes().weight) {
            return false;
        }
        if (connections1[i].get_attributes().enabled != connections2[i].get_attributes().enabled) {
            return false;
        }
        if (connections1[i].get_sourceNodeIndex() != connections2[i].get_sourceNodeIndex()) {
            return false;
        }
        if (connections1[i].get_targetNodeIndex() != connections2[i].get_targetNodeIndex()) {
            return false;
        }
    }
    
    return true;
}

}