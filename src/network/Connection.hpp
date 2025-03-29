// Connection.hpp
#pragma once
#include <memory>
#include <stdexcept>
#include "core/ActivationGene.hpp"
#include "core/ActivationFunction.hpp"

namespace neat {
namespace network {

class Node;
using NodePtr = std::shared_ptr<Node>;

class Connection {
public:
    Connection(NodePtr source, NodePtr target, double weight, bool enabled, 
              core::ActivationGene activation)
        : source(source)
        , target(target)
        , weight(weight)
        , enabled(enabled)
        , activation(activation) {
        if (!source || !target) {
            throw std::invalid_argument("Source and target nodes must be valid");
        }
    }

    // Core accessors
    NodePtr getSource() const { return source; }
    NodePtr getTarget() const { return target; }
    double getWeight() const { return weight; }
    bool isEnabled() const { return enabled; }
    const core::ActivationGene& getActivation() const { return activation; }

    // Network-focused accessors
    int32_t getSourceId() const;// { return source->getId(); }
    int32_t getTargetId() const;// { return target->getId(); }

    // Modifiers
    void setWeight(double w) { weight = w; }
    void setEnabled(bool e) { enabled = e; }
    void setActivation(const core::ActivationGene& act) { activation = act; }

    // Utility methods
    double getRawInput() const;// {
    //     return target->getValue() * weight;
    // }

    double activate(double input) const {
        auto activationFunc = core::ActivationFunction::getFunction(activation.getType());
        return activationFunc(input);
    }

private:
    NodePtr source;
    NodePtr target;
    double weight;
    bool enabled;
    core::ActivationGene activation;
};

}
}
