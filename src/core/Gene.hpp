// gene.hpp
#pragma once
#include <cstdint>
#include "ActivationGene.hpp"

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
    Gene(int32_t input, int32_t output, double w, bool en, int32_t innov)
        : inputNode(input)
        , outputNode(output)
        , weight(w)
        , enabled(en)
        , innovation(innov)
        , activation() {}
        
    bool operator==(const Gene& other) const {
        return innovation == other.innovation;
    }
    
    bool operator<(const Gene& other) const {
        return innovation < other.innovation;
    }
    // In Gene.hpp, add:
    static double getActivationCompatibility(const Gene& gene1, const Gene& gene2) {
        // Return 0.0 if same activation, 1.0 if different
        return gene1.activation.getType() == gene2.activation.getType() ? 0.0 : 1.0;
}

};

}
}
