// network.hpp
#pragma once
#include <map>
#include <stdexcept>
#include "Node.hpp"
#include "Activation.hpp"
#include "core/Gene.hpp"

namespace neat {
namespace network {

class Network {
public:
    struct Config {
        activation::EActivationFunction activationFunction = activation::EActivationFunction::SIGMOID;
        bool allowRecurrent = false;
        double biasValue = 1.0;
    };
    
    Network(const Config& config = Config{})
        : config(config)
        , activationFn(activation::ActivationFunctions::getFunction(config.activationFunction)) {}
    
    void addNode(int32_t id, core::ENodeType type) {
        auto node = std::make_shared<Node>(id, type, activationFn);
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
    
    void addConnection(int32_t fromId, int32_t toId, double weight, bool enabled = true) {
        auto from = fromId == -1 ? biasNode : nodes[fromId];
        auto to = nodes[toId];
        
        if (!from || !to) {
            throw std::runtime_error("Invalid node IDs in addConnection");
        }
        
        if (!config.allowRecurrent && wouldCreateCycle(fromId, toId)) {
            throw std::runtime_error("Attempting to create cycle in feed-forward network");
        }
        
        from->addOutput(to, weight, enabled);
        to->addInput(from, weight, enabled);
        
        if (enabled) {
            connections.emplace_back(fromId, toId, weight);
        }
    }
    
    std::vector<double> activate(const std::vector<double>& inputs) {
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
    
    bool validate() const {
        try {
            // Verify node consistency
            for (const auto& [id, node] : nodes) {
                if (id != node->getId()) {
                    return false;
                }
            }
            
            // Verify connection consistency
            for (const auto& [fromId, toId, weight] : connections) {
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

private:
    struct Connection {
        int32_t fromId;
        int32_t toId;
        double weight;
        
        Connection(int32_t f, int32_t t, double w)
            : fromId(f), toId(t), weight(w) {}
    };
    
    Config config;
    ActivationFunction activationFn;
    std::map<int32_t, NodePtr> nodes;
    std::vector<NodePtr> inputNodes;
    std::vector<NodePtr> hiddenNodes;
    std::vector<NodePtr> outputNodes;
    NodePtr biasNode;
    std::vector<Connection> connections;
    
    std::vector<int32_t> getTopologicalOrder() const;
    bool hasCycles() const;
    bool wouldCreateCycle(int32_t fromId, int32_t toId) const;
};

} // namespace network
} // namespace neat
