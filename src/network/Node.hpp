// node.hpp
#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <string>
#include <functional>
#include "core/Gene.hpp"
#include "Connection.hpp"
#include "core/ActivationGene.hpp"
#include "core/ActivationFunction.hpp"

namespace neat {
namespace network {


class Node {
public:
    Node(int32_t id, core::ENodeType type)
        : id(id)
        , type(type)
        , value(0.0) {}

    // Core accessors
    int32_t getId() const { return id; }
    core::ENodeType getType() const { return type; }
    double getValue() const { return value; }
    void setValue(double v) { value = v; }

    // Connection management
    void addConnection(const Connection& conn) {
        if (conn.getSource()->getId() == id) {
            outgoingConnections.push_back(conn);
        } else if (conn.getTarget()->getId() == id) {
            incomingConnections.push_back(conn);
        } else {
            throw std::invalid_argument("Connection does not involve this node");
        }
    }

    void removeConnection(int32_t otherId) {
        auto removeIf = [otherId](const Connection& conn) {
            return conn.getSourceId() == otherId || conn.getTargetId() == otherId;
        };
        outgoingConnections.erase(
            std::remove_if(outgoingConnections.begin(), outgoingConnections.end(), removeIf),
            outgoingConnections.end()
        );
        incomingConnections.erase(
            std::remove_if(incomingConnections.begin(), incomingConnections.end(), removeIf),
            incomingConnections.end()
        );
    }

    const std::vector<Connection>& getIncoming() const { return incomingConnections; }
    const std::vector<Connection>& getOutgoing() const { return outgoingConnections; }

    bool isInput() const {
        return type == core::ENodeType::INPUT || type == core::ENodeType::BIAS;
    }

private:
    int32_t id;
    core::ENodeType type;
    double value;
    std::vector<Connection> incomingConnections;
    std::vector<Connection> outgoingConnections;
};

}
}
