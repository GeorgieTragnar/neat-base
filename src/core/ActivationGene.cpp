

#include "ActivationGene.hpp"


#include <cstdint>
#include <functional>
#include <random>
#include <vector>
#include <memory>


namespace neat {
namespace core {

ActivationGene ActivationGene::createRandom(const Config& config) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> dis(0.0, 1.0);
	
	double prob = dis(gen);
	EActivationType type;
	
	if (prob < config.basicTierProb) {
		// Basic tier
		type = static_cast<EActivationType>(gen() % 3);
	} else if (prob < config.basicTierProb + config.advancedTierProb) {
		// Advanced tier
		type = static_cast<EActivationType>(3 + (gen() % 2));
	} else {
		// Experimental tier
		type = static_cast<EActivationType>(5 + (gen() % 2));
	}
	
	return ActivationGene(type);
}

void ActivationGene::mutate(const Config& config) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> dis(0.0, 1.0);
	
	if (dis(gen) < config.mutationRate) {
		int currentTier = getTier(type);
		int newTier = currentTier;
		
		// Decide whether to change tiers
		if (dis(gen) < config.interTierMutationRate) {
			// Move up or down one tier, wrapping around
			newTier = (currentTier + (dis(gen) < 0.5 ? -1 : 1)) % 3;
			if (newTier < 0) newTier = 2;
		}
		
		// Select new activation within tier
		switch (newTier) {
			case 0: // Basic
				type = static_cast<EActivationType>(gen() % 3);
				break;
			case 1: // Advanced
				type = static_cast<EActivationType>(3 + (gen() % 2));
				break;
			case 2: // Experimental
				type = static_cast<EActivationType>(5 + (gen() % 2));
				break;
		}
	}
}

}
}