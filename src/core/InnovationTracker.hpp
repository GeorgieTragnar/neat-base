// innovation.hpp
#pragma once
#include <cstdint>
#include <vector>
#include <unordered_map>

namespace neat {
namespace core {

class InnovationTracker {
public:
    // Core NEAT paper method: Get next global innovation number
    static int32_t getNextInnovation() noexcept {
        return nextInnovation++;
    }
    
    // Reset innovation tracking (useful between evolutionary runs)
    static void reset() {
        nextInnovation = 0;
        currentGenerationInnovations.clear();
    }

    // Track innovations within current generation to ensure consistent numbering
    static int32_t trackGenerationInnovation(
        int32_t inputNode, 
        int32_t outputNode
    ) {
        // Generate unique key for connection
        auto key = std::make_pair(inputNode, outputNode);
        
        // Check if this exact connection was already innovated this generation
        auto it = currentGenerationInnovations.find(key);
        if (it != currentGenerationInnovations.end()) {
            return it->second;
        }
        
        // If not, assign new innovation number and track it
        int32_t innovationNumber = getNextInnovation();
        currentGenerationInnovations[key] = innovationNumber;
        
        return innovationNumber;
    }

    
private:
    // Global innovation number tracker
    static int32_t nextInnovation;

    // Custom hash for pair to use in unordered_map
    struct PairHash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ h2;
        }
    };

    // Track innovations within current generation
    // Ensures same structural mutation gets same innovation number
    static std::unordered_map<
        std::pair<int32_t, int32_t>,  // input and output node
        int32_t,  // innovation number
        PairHash
    > currentGenerationInnovations;
};

}
}
