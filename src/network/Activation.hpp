// activation.hpp
#pragma once
#include <cmath>
#include "Node.hpp"

namespace neat {
namespace network {
namespace activation {

enum class EActivationFunction {
	SIGMOID,
	TANH,
	RELU,
	LEAKY_RELU,
	SOFTPLUS
};

class ActivationFunctions {
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
    
    static ActivationFunction getFunction(const EActivationFunction& funcType) {
		switch (funcType)
		{
			case EActivationFunction::SIGMOID: return sigmoid;
			case EActivationFunction::TANH : return tanh;
			case EActivationFunction::RELU : return relu;
			case EActivationFunction::LEAKY_RELU : return leakyRelu;
			case EActivationFunction::SOFTPLUS : return softplus;
			default: return sigmoid;	
		}
    }
};

} // namespace activation
}
}
