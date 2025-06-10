#include "CompatibilityDistance.hpp"
#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace Operator {

CompatibilityDistanceParams::CompatibilityDistanceParams(float c1, float c2, float c3, float threshold)
    : _c1(c1), _c2(c2), _c3(c3), _threshold(threshold) {
}

uint32_t compatibilityDistance(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const CompatibilityDistanceParams& params) {
    // Iterate through existing species representatives
    for (const auto& [speciesId, representative] : historyTracker->_speciesRepresentatives) {
        const auto& connections1 = genome.get_connectionGenes();
        const auto& connections2 = representative.get_connectionGenes();
        
        // Create maps for efficient lookup by history ID
        std::unordered_map<uint32_t, size_t> conn1_map;
        std::unordered_map<uint32_t, size_t> conn2_map;
        
        for (size_t i = 0; i < connections1.size(); ++i) {
            conn1_map[connections1[i].get_historyID()] = i;
        }
        
        for (size_t i = 0; i < connections2.size(); ++i) {
            conn2_map[connections2[i].get_historyID()] = i;
        }
        
        // Find max history ID in both genomes
        uint32_t max_id1 = connections1.empty() ? 0 : 
            std::max_element(connections1.begin(), connections1.end(),
                [](const ConnectionGene& a, const ConnectionGene& b) {
                    return a.get_historyID() < b.get_historyID();
                })->get_historyID();
        
        uint32_t max_id2 = connections2.empty() ? 0 :
            std::max_element(connections2.begin(), connections2.end(),
                [](const ConnectionGene& a, const ConnectionGene& b) {
                    return a.get_historyID() < b.get_historyID();
                })->get_historyID();
        
        uint32_t max_id = std::max(max_id1, max_id2);
        
        size_t excess = 0;
        size_t disjoint = 0;
        size_t matching = 0;
        float weight_diff_sum = 0.0f;
        
        // Count excess, disjoint, and matching genes
        for (const auto& conn1 : connections1) {
            uint32_t id = conn1.get_historyID();
            auto it2 = conn2_map.find(id);
            
            if (it2 != conn2_map.end()) {
                // Matching gene - calculate weight difference
                matching++;
                weight_diff_sum += std::abs(conn1.get_attributes().weight - 
                                          connections2[it2->second].get_attributes().weight);
            } else {
                // Gene exists in genome1 but not genome2
                if (id > max_id2) {
                    excess++;
                } else {
                    disjoint++;
                }
            }
        }
        
        // Count disjoint/excess genes that exist in genome2 but not genome1
        for (const auto& conn2 : connections2) {
            uint32_t id = conn2.get_historyID();
            if (conn1_map.find(id) == conn1_map.end()) {
                if (id > max_id1) {
                    excess++;
                } else {
                    disjoint++;
                }
            }
        }
        
        // Calculate compatibility distance
        size_t N = std::max(connections1.size(), connections2.size());
        if (N < 20) N = 1; // Small genomes don't normalize
        
        float average_weight_diff = matching > 0 ? weight_diff_sum / matching : 0.0f;
        
        float distance = (params._c1 * excess) / N + (params._c2 * disjoint) / N + params._c3 * average_weight_diff;
        
        if (distance < params._threshold) {
            return speciesId;
        }
    }
    
    // No compatible species found - create new species
    uint32_t newSpeciesId = historyTracker->_nextSpeciesID++;
    historyTracker->_speciesRepresentatives.emplace(newSpeciesId, genome);
    
    return newSpeciesId;
}

}