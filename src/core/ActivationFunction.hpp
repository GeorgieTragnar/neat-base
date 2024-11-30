// ActivationFunction.hpp
#pragma once
#include <cmath>

namespace neat {
namespace core {

class ActivationFunction {
public:
    static double sigmoid(double x) {
        return 1.0 / (1.0 + std::exp(-4.9 * x));
    }
    
    static double tanh(double x) {
        return std::tanh(x);
    }
    
    static double relu(double x) {
        return std::max(0.0, x);
    }
    
    static double leakyRelu(double x) {
        return x > 0 ? x : 0.01 * x;
    }
    
    static double softplus(double x) {
        return std::log1p(std::exp(x));
    }
    
    static double gaussian(double x) {
        return std::exp(-x * x);
    }
    
    static double sine(double x) {
        return std::sin(x);
    }
    
    static std::function<double(double)> getFunction(EActivationType type) {
        switch (type) {
            case EActivationType::SIGMOID: return sigmoid;
            case EActivationType::TANH: return tanh;
            case EActivationType::RELU: return relu;
            case EActivationType::LEAKY_RELU: return leakyRelu;
            case EActivationType::SOFTPLUS: return softplus;
            case EActivationType::GAUSSIAN: return gaussian;
            case EActivationType::SINE: return sine;
            default: return sigmoid;
        }
    }
};

} // namespace core
} // namespace neat
