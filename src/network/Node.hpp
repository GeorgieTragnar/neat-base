// node.hpp
#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <string>
#include <functional>
#include "core/Gene.hpp"

namespace neat {
namespace network {

class Node;
using NodePtr = std::shared_ptr<Node>;
using ActivationFunction = std::function<double(double)>;

class Node {
public:
    struct Connection {
        NodePtr target;
        double weight;
        bool enabled;
        
        Connection(NodePtr t, double w, bool e = true)
            : target(t), weight(w), enabled(e) {}
    };
    
    Node(int32_t id, core::ENodeType type, ActivationFunction actFunc)
        : nodeId(id)
        , ENodeType(type)
        , activation(std::move(actFunc))
        , value(0.0)
        , sum(0.0) {}
    
    void addInput(NodePtr node, double weight, bool enabled = true) {
        inputs.emplace_back(node, weight, enabled);
    }
    
    void addOutput(NodePtr node, double weight, bool enabled = true) {
        outputs.emplace_back(node, weight, enabled);
    }
    
    double activate() {
        if (ENodeType == core::ENodeType::INPUT || ENodeType == core::ENodeType::BIAS) {
            return value;
        }
        
        sum = 0.0;
        for (const auto& conn : inputs) {
            if (conn.enabled) {
                sum += conn.target->getValue() * conn.weight;
            }
        }
        
        value = activation(sum);
        return value;
    }
    
    // Getters and setters
    void setValue(double v) noexcept { value = v; }
    double getValue() const noexcept { return value; }
    int32_t getId() const noexcept { return nodeId; }
    core::ENodeType getType() const noexcept { return ENodeType; }
    
    const std::vector<Connection>& getInputs() const noexcept { return inputs; }
    const std::vector<Connection>& getOutputs() const noexcept { return outputs; }

private:
    int32_t nodeId;
    core::ENodeType ENodeType;
    ActivationFunction activation;
    std::vector<Connection> inputs;
    std::vector<Connection> outputs;
    double value;
    double sum;
};

}
}
