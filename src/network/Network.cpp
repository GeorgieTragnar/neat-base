
#include "Network.hpp"
#include <queue>
#include <set>

namespace neat {
namespace network {

void Network::addNode(int32_t id, core::ENodeType type) {
    auto node = std::make_shared<Node>(id, type);
    nodes[id] = node;
    
    switch (type) {
        case core::ENodeType::INPUT:
            if (inputNodes.size() >= config.inputSize) {
                throw std::runtime_error("Exceeding maximum input size");
            }
            inputNodes.push_back(node);
            break;
        case core::ENodeType::OUTPUT:
            if (outputNodes.size() >= config.outputSize) {
                throw std::runtime_error("Exceeding maximum output size");
            }
            outputNodes.push_back(node);
            break;
        case core::ENodeType::BIAS:
            if (biasNode) {
                throw std::runtime_error("Multiple bias nodes not allowed");
            }
            biasNode = node;
            biasNode->setValue(config.biasValue);
            break;
        default:
            hiddenNodes.push_back(node);
            break;
    }
}

void Network::addConnection(int32_t fromId, int32_t toId, double weight, bool enabled, core::EActivationType actType) {
	auto sourceNode = fromId == -1 ? biasNode : nodes[fromId];
	auto targetNode = nodes[toId];
	
	if (!sourceNode || !targetNode) {
		throw std::runtime_error("Invalid node IDs in addConnection");
	}
	
	if (!config.allowRecurrent && wouldCreateCycle(fromId, toId)) {
		throw std::runtime_error("Attempting to create cycle in feed-forward network");
	}
	
    Connection conn(sourceNode, targetNode, weight, enabled, core::ActivationGene(actType));

    sourceNode->addConnection(conn);
    targetNode->addConnection(conn);
}

std::vector<double> Network::activate(const std::vector<double>& inputs) {
    if (inputs.size() != config.inputSize) {
        throw std::runtime_error("Invalid input size: " + 
            std::to_string(inputs.size()) + 
            " expected: " + std::to_string(config.inputSize));
    }
	
	// Set input values
	for (size_t i = 0; i < inputs.size(); ++i) {
		inputNodes[i]->setValue(inputs[i]);
	}
	
	// Activate nodes in topological order
	auto order = getTopologicalOrder();
	for (const auto& nodeId : order) {
		activateNode(nodes[nodeId]);
	}
	
	// Collect outputs
	std::vector<double> outputs;
	outputs.reserve(outputNodes.size());
	for (const auto& node : outputNodes) {
		outputs.push_back(node->getValue());
	}
	
	return outputs;
}

void Network::activateNode(NodePtr node) {
    // Skip activation for input nodes (including bias)
    if (node->isInput()) {
        return;
    }
    
    double sum = 0.0;
    for (const auto& conn : node->getIncoming()) {
        if (conn.isEnabled()) {
            double input = conn.getRawInput();
            sum += conn.activate(input);
        }
    }
    
    node->setValue(sum);
}

double Network::activateConnection(const Connection& conn) {
    double rawInput = conn.getRawInput();
    auto activationFunc = core::ActivationFunction::getFunction(conn.getActivation().getType());
    return activationFunc(rawInput);
}

bool Network::validate() const {
	try {
		// Verify node consistency
		for (const auto& [id, node] : nodes) {
			if (id != node->getId()) {
				return false;
			}
		}
        
        // Get all connections for validation
        auto allConnections = getAllConnections();
        
        // Verify connection consistency
        for (const auto& conn : allConnections) {
            auto from = conn.get().getSource();
            auto to = conn.get().getTarget();
            
            // Verify nodes exist in our map
            if (from->getId() != -1 && nodes.find(from->getId()) == nodes.end()) return false;
            if (nodes.find(to->getId()) == nodes.end()) return false;
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

// Helper method to get all connections in the network
std::vector<std::reference_wrapper<Connection>> Network::getAllConnections() const {
    std::vector<std::reference_wrapper<Connection>> allConnections;
    
    // Collect all connections from nodes
    for (const auto& [id, node] : nodes) {
        // Add outgoing connections
        for (const auto& conn : node->getOutgoing()) {
            // Only add each connection once (when it's outgoing from a node)
            allConnections.push_back(std::ref(const_cast<Connection&>(conn)));
        }
    }
    
    // Also check bias node if it exists
    if (biasNode) {
        for (const auto& conn : biasNode->getOutgoing()) {
            allConnections.push_back(std::ref(const_cast<Connection&>(conn)));
        }
    }
    
    return allConnections;
}

// Implement topological sorting using Kahn's algorithm
std::vector<int32_t> Network::getTopologicalOrder() const {
    std::vector<int32_t> result;
    std::map<int32_t, int32_t> inDegree;
    std::queue<int32_t> queue;
    
    // Initialize in-degree for all nodes
    for (const auto& [id, _] : nodes) {
        inDegree[id] = 0;
    }
    
    // Calculate in-degree for each node using actual connections
    auto allConnections = getAllConnections();
    for (const auto& connRef : allConnections) {
        const auto& conn = connRef.get();
        if (conn.isEnabled()) {
            int32_t targetId = conn.getTargetId();
            int32_t sourceId = conn.getSourceId();
            
            // Only count connections that aren't from bias nodes
            // when calculating in-degree
            if (sourceId != -1) {  // -1 is the bias node ID
                inDegree[targetId]++;
            }
        }
    }
    
    // Add bias node to result first (if it exists)
    if (biasNode) {
        result.push_back(-1);  // Bias node ID is -1
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
        
        // Get the current node
        auto nodeIt = nodes.find(current);
        if (nodeIt != nodes.end()) {
            // Process outgoing connections
            for (const auto& conn : nodeIt->second->getOutgoing()) {
                if (conn.isEnabled()) {
                    int32_t targetId = conn.getTargetId();
                    inDegree[targetId]--;
                    if (inDegree[targetId] == 0) {
                        queue.push(targetId);
                    }
                }
            }
        }
    }
    
    // If result size doesn't match node count (+1 for bias if present), there's a cycle
    size_t expectedSize = nodes.size() + (biasNode ? 1 : 0);
    return result.size() == expectedSize ? result : std::vector<int32_t>();
}

// Check for cycles using DFS
bool Network::hasCycles() const {
    std::set<int32_t> visited;
    std::set<int32_t> recursionStack;
    
    // Helper function for DFS
    std::function<bool(int32_t)> hasCycleUtil = [&](int32_t nodeId) -> bool {
        visited.insert(nodeId);
        recursionStack.insert(nodeId);
        
        // Get the node
        auto nodeIt = nodes.find(nodeId);
        if (nodeIt != nodes.end()) {
            // Check all outgoing connections
            for (const auto& conn : nodeIt->second->getOutgoing()) {
                if (conn.isEnabled()) {
                    int32_t neighbor = conn.getTargetId();
                    
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
        
        // Get the start node
        auto nodeIt = nodes.find(start);
        if (nodeIt != nodes.end()) {
            // Check all outgoing connections
            for (const auto& conn : nodeIt->second->getOutgoing()) {
                if (conn.isEnabled()) {
                    int32_t nextId = conn.getTargetId();
                    if (visited.find(nextId) == visited.end()) {
                        if (hasPath(nextId, target)) {
                            return true;
                        }
                    }
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
