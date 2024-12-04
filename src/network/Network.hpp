// network.hpp
#pragma once
#include <map>
#include <stdexcept>
#include "Node.hpp"
#include "core/Gene.hpp"

namespace neat {
namespace network {

class Network {
public:
    struct Config {
        // Core network parameters
        int32_t inputSize;
        int32_t outputSize;
        bool allowRecurrent;
        double biasValue;
        core::ActivationGene::Config activationConfig;

        Config() = delete;
        Config(int32_t inputs, 
               int32_t outputs,
               const core::ActivationGene::Config& actConfig,
               bool recurrent = false,
               double bias = 1.0)
            : inputSize(inputs)
            , outputSize(outputs)
            , allowRecurrent(recurrent)
            , biasValue(bias)
            , activationConfig(actConfig) {}
    };

    explicit Network(const Config& config)
        : config(config) {
        // Pre-allocate vectors for performance
        inputNodes.reserve(config.inputSize);
        outputNodes.reserve(config.outputSize);
        hiddenNodes.reserve(10);
    }
    
    void addNode(int32_t id, core::ENodeType type, core::EActivationType actType = core::EActivationType::SIGMOID);
    void addConnection(int32_t fromId, int32_t toId, double weight, bool enabled = true, core::EActivationType actType = core::EActivationType::SIGMOID);
    std::vector<double> activate(const std::vector<double>& inputs);
    bool validate() const;

    const Config getConfig() const { return config; }

private:
    struct Connection {
        int32_t fromId;
        int32_t toId;
        double weight;
        core::EActivationType actType;
        
        Connection(int32_t f, int32_t t, double w, core::EActivationType type)
            : fromId(f), toId(t), weight(w), actType(type) {}
    };
    
    Config config;
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
