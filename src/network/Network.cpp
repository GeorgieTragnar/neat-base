
#include "Network.hpp"
#include <queue>
#include <set>

namespace neat {
namespace network {

void Network::addNode(int32_t id, core::ENodeType type, core::EActivationType actType) {
	auto node = std::make_shared<Node>(id, type, actType);
	nodes[id] = node;
	
	switch (type) {
		case core::ENodeType::INPUT:
			inputNodes.push_back(node);
			break;
		case core::ENodeType::OUTPUT:
			outputNodes.push_back(node);
			break;
		case core::ENodeType::BIAS:
			biasNode = node;
			biasNode->setValue(config.biasValue);
			break;
		default:
			hiddenNodes.push_back(node);
			break;
	}
}

void Network::addConnection(int32_t fromId, int32_t toId, double weight, bool enabled, core::EActivationType actType) {
	auto from = fromId == -1 ? biasNode : nodes[fromId];
	auto to = nodes[toId];
	
	if (!from || !to) {
		throw std::runtime_error("Invalid node IDs in addConnection");
	}
	
	if (!config.allowRecurrent && wouldCreateCycle(fromId, toId)) {
		throw std::runtime_error("Attempting to create cycle in feed-forward network");
	}
	
	from->addOutput(to, weight, enabled, actType);
	to->addInput(from, weight, enabled, actType);
	
	if (enabled) {
		connections.emplace_back(fromId, toId, weight, actType);
	}
}

std::vector<double> Network::activate(const std::vector<double>& inputs) {
	if (inputs.size() != inputNodes.size()) {
		throw std::runtime_error("Invalid input size");
	}
	
	// Set input values
	for (size_t i = 0; i < inputs.size(); ++i) {
		inputNodes[i]->setValue(inputs[i]);
	}
	
	// Activate nodes in topological order
	auto order = getTopologicalOrder();
	for (const auto& nodeId : order) {
		nodes[nodeId]->activate();
	}
	
	// Collect outputs
	std::vector<double> outputs;
	outputs.reserve(outputNodes.size());
	for (const auto& node : outputNodes) {
		outputs.push_back(node->getValue());
	}
	
	return outputs;
}

bool Network::validate() const {
	try {
		// Verify node consistency
		for (const auto& [id, node] : nodes) {
			if (id != node->getId()) {
				return false;
			}
		}
		
		// Verify connection consistency
		for (const auto& [fromId, toId, weight, actType] : connections) {
			auto from = fromId == -1 ? biasNode : nodes.at(fromId);
			auto to = nodes.at(toId);
			
			bool connectionFound = false;
			for (const auto& conn : from->getOutputs()) {
				if (conn.target == to && conn.weight == weight) {
					connectionFound = true;
					break;
				}
			}
			
			if (!connectionFound) return false;
		}
		
		// Verify no cycles (if feed-forward)
		if (!config.allowRecurrent && hasCycles()) {
			return false;
		}
		
		return true;
	} catch (const std::exception&) {
		return false;
	}
}

// Implement topological sorting using Kahn's algorithm
std::vector<int32_t> neat::network::Network::getTopologicalOrder() const {
    std::vector<int32_t> result;
    std::map<int32_t, int32_t> inDegree;
    std::queue<int32_t> queue;
    
    // Initialize in-degree for all nodes
    for (const auto& [id, _] : nodes) {
        inDegree[id] = 0;
    }
    
    // Calculate in-degree for each node
    for (const auto& conn : connections) {
        if (conn.toId != -1) { // Skip bias node
            inDegree[conn.toId]++;
        }
    }
    
    // Add nodes with no incoming edges to queue
    for (const auto& [id, _] : nodes) {
        if (inDegree[id] == 0) {
            queue.push(id);
        }
    }
    
    // Process nodes
    while (!queue.empty()) {
        int32_t current = queue.front();
        queue.pop();
        result.push_back(current);
        
        // Decrease in-degree for all neighbors
        for (const auto& conn : connections) {
            if (conn.fromId == current) {
                inDegree[conn.toId]--;
                if (inDegree[conn.toId] == 0) {
                    queue.push(conn.toId);
                }
            }
        }
    }
    
    // If result size doesn't match node count, there's a cycle
    return result.size() == nodes.size() ? result : std::vector<int32_t>();
}

// Check for cycles using DFS
bool Network::hasCycles() const {
    std::set<int32_t> visited;
    std::set<int32_t> recursionStack;
    
    // Helper function for DFS
    std::function<bool(int32_t)> hasCycleUtil = [&](int32_t nodeId) -> bool {
        visited.insert(nodeId);
        recursionStack.insert(nodeId);
        
        // Check all outgoing connections
        for (const auto& conn : connections) {
            if (conn.fromId == nodeId) {
                int32_t neighbor = conn.toId;
                
                if (visited.find(neighbor) == visited.end()) {
                    if (hasCycleUtil(neighbor)) {
                        return true;
                    }
                }
                else if (recursionStack.find(neighbor) != recursionStack.end()) {
                    return true;
                }
            }
        }
        
        recursionStack.erase(nodeId);
        return false;
    };
    
    // Check from each unvisited node
    for (const auto& [id, _] : nodes) {
        if (visited.find(id) == visited.end()) {
            if (hasCycleUtil(id)) {
                return true;
            }
        }
    }
    
    return false;
}

// Check if adding a new connection would create a cycle
bool Network::wouldCreateCycle(int32_t fromId, int32_t toId) const {
    if (fromId == toId) {
        return true;  // Self-connection is a cycle
    }
    
    // Use DFS to check if there's already a path from toId to fromId
    // If such a path exists, adding fromId->toId would create a cycle
    std::set<int32_t> visited;
    
    std::function<bool(int32_t, int32_t)> hasPath = [&](int32_t start, int32_t target) -> bool {
        if (start == target) {
            return true;
        }
        
        visited.insert(start);
        
        // Check all outgoing connections from current node
        for (const auto& conn : connections) {
            if (conn.fromId == start && visited.find(conn.toId) == visited.end()) {
                if (hasPath(conn.toId, target)) {
                    return true;
                }
            }
        }
        
        return false;
    };
    
    // If there's a path from toId to fromId, adding fromId->toId would create a cycle
    return hasPath(toId, fromId);
}

}
}
