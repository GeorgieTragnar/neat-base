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
    
    Network(const Config& config)
        : config(config)
        , activationFn(activation::ActivationFunctions::getFunction(config.activationFunction)) {}
    
    void addNode(int32_t id, core::ENodeType type);
    void addConnection(int32_t fromId, int32_t toId, double weight, bool enabled = true);
    std::vector<double> activate(const std::vector<double>& inputs);
    bool validate() const;

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
