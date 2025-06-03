#pragma once

namespace Analysis {

class FitnessResult {
public:
    FitnessResult(double value) : _value(value) {}
    
    bool isBetterThan(const FitnessResult& other) const {
        return _value > other._value;
    }
    
    bool isEqualTo(const FitnessResult& other) const {
        return _value == other._value;
    }
    
    bool operator>(const FitnessResult& other) const {
        return isBetterThan(other);
    }
    
    bool operator==(const FitnessResult& other) const {
        return isEqualTo(other);
    }

private:
    double _value;
};

}