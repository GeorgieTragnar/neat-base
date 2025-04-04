// gene.hpp
#pragma once
#include <cstdint>
#include <limits>

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
    
    Gene() = delete;
    Gene(int32_t input, int32_t output, double w, bool en, int32_t innov)
        : inputNode(input)
        , outputNode(output)
        , weight(w)
        , enabled(en)
        , innovation(innov) {}
        
    bool operator==(const Gene& other) const {
        return innovation == other.innovation;
    }
    
    bool operator<(const Gene& other) const {
        return innovation < other.innovation;
    }
};

}
}
