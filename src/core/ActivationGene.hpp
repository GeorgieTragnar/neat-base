// ActivationGene.hpp
#pragma once

#include <stdexcept>

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
			// Mutation parameters
			double mutationRate;  // Chance of mutating activation
			double interTierMutationRate; // Chance of changing tiers when mutating
			
			// Tier probabilities
			double basicTierProb;     // For initial random assignment
			double advancedTierProb;  // For initial random assignment
			// Experimental tier = 1.0 - (basic + advanced)
			
			// Speciation parameter
			double compatibilityWeight; // How much activation differences matter in speciation
			
            Config() = delete;
			Config(double mutRate = 0.1,
				   double interTierRate = 0.2,
				   double basicProb = 0.8,
				   double advancedProb = 0.15,
				   double compatWeight = 1.0) {
				
				if (basicProb + advancedProb > 1.0) {
					throw std::runtime_error("Tier probabilities exceed 1.0");
				}
				
				mutationRate = mutRate;
				interTierMutationRate = interTierRate;
				basicTierProb = basicProb;
				advancedTierProb = advancedProb;
				compatibilityWeight = compatWeight;
			}
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
