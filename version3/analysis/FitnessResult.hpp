#pragma once
#include <memory>

namespace Analysis {

class FitnessResultInterface {
public:
    virtual ~FitnessResultInterface() = default;
    
    virtual bool isBetterThan(const FitnessResultInterface& other) const = 0;
    virtual bool isEqualTo(const FitnessResultInterface& other) const = 0;
    
    virtual std::unique_ptr<FitnessResultInterface> clone() const = 0;
};

}