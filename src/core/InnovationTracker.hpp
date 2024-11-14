// innovation.hpp
#pragma once
#include <cstdint>

namespace neat {
namespace core {

class InnovationTracker {
public:
    static int32_t getNextInnovation() noexcept {
        return nextInnovation++;
    }
    
    static void reset() noexcept {
        nextInnovation = 0;
    }
    
private:
    static int32_t nextInnovation;
};

}
}
