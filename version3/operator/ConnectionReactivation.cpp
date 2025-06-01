#include "ConnectionReactivation.hpp"
#include <random>
#include <cassert>
#include <algorithm>

using namespace Operator;

ConnectionReactivationParams::ConnectionReactivationParams(SelectionStrategy strategy)
    : _selectionStrategy(strategy)
{
    // No parameter validation needed - enum is always valid
}

namespace {
    // Helper to find all disabled connection indices
    std::vector<size_t> findDisabledConnectionIndices(const std::vector<ConnectionGene>& connections) {
        std::vector<size_t> disabledIndices;
        
        for (size_t i = 0; i < connections.size(); ++i) {
            if (!connections[i].get_attributes().enabled) {
                disabledIndices.push_back(i);
            }
        }
        
        return disabledIndices;
    }
    
    // Helper to select connection index based on strategy
    size_t selectConnectionIndex(const std::vector<size_t>& disabledIndices,
                                const std::vector<ConnectionGene>& connections,
                                ConnectionReactivationParams::SelectionStrategy strategy) {
        
        switch (strategy) {
            case ConnectionReactivationParams::SelectionStrategy::RANDOM: {
                static std::random_device rd;
                static std::mt19937 gen(rd());
                std::uniform_int_distribution<size_t> dist(0, disabledIndices.size() - 1);
                return disabledIndices[dist(gen)];
            }
            
            case ConnectionReactivationParams::SelectionStrategy::OLDEST_FIRST: {
                // Find disabled connection with lowest innovation number
                size_t oldestIndex = disabledIndices[0];
                uint32_t lowestInnovation = connections[oldestIndex].get_historyID();
                
                for (size_t idx : disabledIndices) {
                    if (connections[idx].get_historyID() < lowestInnovation) {
                        lowestInnovation = connections[idx].get_historyID();
                        oldestIndex = idx;
                    }
                }
                return oldestIndex;
            }
            
            case ConnectionReactivationParams::SelectionStrategy::NEWEST_FIRST: {
                // Find disabled connection with highest innovation number
                size_t newestIndex = disabledIndices[0];
                uint32_t highestInnovation = connections[newestIndex].get_historyID();
                
                for (size_t idx : disabledIndices) {
                    if (connections[idx].get_historyID() > highestInnovation) {
                        highestInnovation = connections[idx].get_historyID();
                        newestIndex = idx;
                    }
                }
                return newestIndex;
            }
        }
        
        // Should never reach here
        assert(false);
        return 0;
    }
}

Genome Operator::connectionReactivation(const Genome& genome, const ConnectionReactivationParams& params) {
    // Create a copy of the genome
    Genome mutatedGenome = genome;
    
    auto& connections = mutatedGenome.get_connectionGenes();
    
    // Find all disabled connections
    auto disabledIndices = findDisabledConnectionIndices(connections);
    
    // Assert if no disabled connections available (user error - shouldn't call operator in this state)
    assert(!disabledIndices.empty());
    
    // Select connection to reactivate based on strategy
    size_t selectedIndex = selectConnectionIndex(disabledIndices, connections, params._selectionStrategy);
    
    // Reactivate the selected connection
    auto& selectedConnection = connections[selectedIndex];
    auto& mutableAttributes = const_cast<ConnectionGeneAttributes&>(selectedConnection.get_attributes());
    mutableAttributes.enabled = true;
    
    return mutatedGenome;
}