// ActivationGene.hpp
#pragma once

namespace neat {
namespace core {

enum class EActivationType {
    // Basic Tier (80% probability)
    SIGMOID,
    TANH,
    RELU,
    
    // Advanced Tier (15% probability)
    LEAKY_RELU,
    SOFTPLUS,
    
    // Experimental Tier (5% probability)
    GAUSSIAN,
    SINE,
    
    COUNT
};

class ActivationGene {
public:
    struct Config {
        double mutationRate = 0.1;
        double interTierMutationRate = 0.2;
        double basicTierProb = 0.8;
        double advancedTierProb = 0.15;
		double compatibilityWeight = 0.5;
    };
    
    ActivationGene() : type(EActivationType::SIGMOID) {}
    explicit ActivationGene(EActivationType t) : type(t) {}
    
    static ActivationGene createRandom(const Config& config);
    
    void mutate(const Config& config);
    
    static int getTier(EActivationType type) {
        if (type <= EActivationType::RELU) return 0;
        if (type <= EActivationType::SOFTPLUS) return 1;
        return 2;
    }
    
    EActivationType getType() const { return type; }
    
private:
    EActivationType type;
};

} // namespace core
} // namespace neat
