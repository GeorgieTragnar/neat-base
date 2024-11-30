// node.hpp
#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <string>
#include <functional>
#include "core/Gene.hpp"
#include "core/ActivationGene.hpp"
#include "core/ActivationFunction.hpp"

namespace neat {
namespace network {

class Node;
using NodePtr = std::shared_ptr<Node>;

class Node {
public:
    struct Connection {
        NodePtr target;
        double weight;
        bool enabled;
        core::EActivationType actType;
        
        Connection(NodePtr t, double w, bool e = true, core::EActivationType act = core::EActivationType::SIGMOID)
            : target(t), weight(w), enabled(e), actType(act) {}
    };
    
    Node(int32_t id, core::ENodeType type, core::EActivationType actType)
        : nodeId(id)
        , ENodeType(type)
        , activationGene(actType)
        , value(0.0)
        , sum(0.0) {
            activationFunc = core::ActivationFunction::getFunction(activationGene.getType());
        }
    
    void addInput(NodePtr node, double weight, bool enabled = true,
                core::EActivationType actType = core::EActivationType::SIGMOID) {
        inputs.emplace_back(node, weight, enabled, actType);
    }
    
    void addOutput(NodePtr node, double weight, bool enabled = true,
                core::EActivationType actType = core::EActivationType::SIGMOID) {
        outputs.emplace_back(node, weight, enabled, actType);
    }
    
    double activate() {
        if (ENodeType == core::ENodeType::INPUT || ENodeType == core::ENodeType::BIAS) {
            return value;
        }
        
        sum = 0.0;
        for (const auto& conn : inputs) {
            if (conn.enabled) {
                // Use the connection's activation function 
                auto inputValue = conn.target->getValue();
                auto actFunc = core::ActivationFunction::getFunction(conn.actType);
                sum += actFunc(inputValue * conn.weight);
            }
        }

        return value = sum;
    }

    void mutateActivation(const core::ActivationGene::Config& config) {
        activationGene.mutate(config);
        activationFunc = core::ActivationFunction::getFunction(activationGene.getType());
    }
    
    // Getters and setters
    void setValue(double v) noexcept { value = v; }
    double getValue() const noexcept { return value; }
    int32_t getId() const noexcept { return nodeId; }
    core::ENodeType getType() const noexcept { return ENodeType; }
    const core::ActivationGene& getActivationGene() const noexcept { return activationGene; }
    
    const std::vector<Connection>& getInputs() const noexcept { return inputs; }
    const std::vector<Connection>& getOutputs() const noexcept { return outputs; }

private:
    int32_t nodeId;
    core::ENodeType ENodeType;
    core::ActivationGene activationGene;
    std::function<double(double)> activationFunc;
    std::vector<Connection> inputs;
    std::vector<Connection> outputs;
    double value;
    double sum;
};

}
}
