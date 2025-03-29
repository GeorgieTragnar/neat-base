// gene.hpp
#pragma once
#include <cstdint>
#include <map>
#include "ActivationGene.hpp"
#include "network/Connection.hpp"

namespace neat {
namespace core {

enum class ENodeType {
    INPUT,
    HIDDEN,
    OUTPUT,
    BIAS
};

struct Gene {
    int32_t inputNode;
    int32_t outputNode;
    double weight;
    bool enabled;
    int32_t innovation;
    core::ActivationGene activation;
    
    Gene() = default;
    Gene(int32_t input, int32_t output, double w, bool en, int32_t innov, core::ActivationGene act)
        : inputNode(input)
        , outputNode(output)
        , weight(w)
        , enabled(en)
        , innovation(innov)
        , activation(act) {}
        
    bool operator==(const Gene& other) const {
        return innovation == other.innovation;
    }
    
    bool operator<(const Gene& other) const {
        return innovation < other.innovation;
    }
    
    static double getActivationCompatibility(const Gene& gene1, const Gene& gene2) {
        // Return 0.0 if same activation, 1.0 if different
        return gene1.activation.getType() == gene2.activation.getType() ? 0.0 : 1.0;
    }

    network::Connection toConnection(
        const std::map<int32_t, network::NodePtr>& nodeMap) const {
        auto source = nodeMap.at(inputNode);
        auto target = nodeMap.at(outputNode);
        return network::Connection(source, target, weight, enabled, activation);
    }

    static Gene fromConnection(const network::Connection& conn, int32_t innovation) {
        return Gene{
            conn.getSourceId(),
            conn.getTargetId(),
            conn.getWeight(),
            conn.isEnabled(),
            innovation,
            conn.getActivation()
        };
    }

};

}
}
