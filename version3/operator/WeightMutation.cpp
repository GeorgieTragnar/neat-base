#include "WeightMutation.hpp"
// #include "../logger/Logger.hpp"
#include <random>
#include <cassert>
#include <algorithm>

using namespace Operator;

WeightMutationParams::WeightMutationParams(double perturbationRate,
                                          double replacementRate, 
                                          double perturbationStrength,
                                          double weightRange,
                                          MutationType mutationType)
    : _perturbationRate(perturbationRate)
    , _replacementRate(replacementRate)
    , _perturbationStrength(perturbationStrength)
    , _weightRange(weightRange)
    , _mutationType(mutationType)
{
    // Parameter validation
    assert(perturbationRate >= 0.0 && perturbationRate <= 1.0);
    assert(replacementRate >= 0.0 && replacementRate <= 1.0);
    assert(perturbationStrength > 0.0);
    assert(weightRange > 0.0);
}

Genome Operator::weightMutation(const Genome& genome, const WeightMutationParams& params) {
    // Create a copy of the genome
    Genome mutatedGenome = genome;
    
    // Assert that deltas are empty before operation
    assert(mutatedGenome.get_connectionGeneDeltas().empty() && "Connection deltas must be empty before weight mutation");
    assert(mutatedGenome.get_nodeGeneDeltas().empty() && "Node deltas must be empty before weight mutation");
    
    // Get mutable reference to connection genes and deltas
    auto& connectionGenes = mutatedGenome.get_connectionGenes();
    auto& connectionDeltas = mutatedGenome.get_connectionGeneDeltas();
    
    if (connectionGenes.empty()) {
        // Debug: Log when no connections exist
        // static auto logger = LOGGER("operator.WeightMutation");
        // LOG_DEBUG("WeightMutation: No connections to mutate, returning unchanged genome");
        return mutatedGenome; // No connections to mutate
    }
    
    // Setup random number generators
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    std::normal_distribution<double> gaussianDist(0.0, params._perturbationStrength);
    std::uniform_real_distribution<double> weightDist(-params._weightRange, params._weightRange);
    
    // Mutate each connection
    for (auto& connection : connectionGenes) {
        // Only mutate enabled connections
        if (!connection.get_attributes().enabled) {
            continue;
        }
        
        auto attributes = connection.get_attributes();
        bool shouldMutate = false;
        
        // Determine mutation type based on strategy
        switch (params._mutationType) {
            case WeightMutationParams::MutationType::PERTURBATION_ONLY:
                if (uniformDist(gen) < params._perturbationRate) {
                    // Perturb weight with Gaussian noise
                    attributes.weight += static_cast<float>(gaussianDist(gen));
                    shouldMutate = true;
                }
                break;
                
            case WeightMutationParams::MutationType::REPLACEMENT_ONLY:
                if (uniformDist(gen) < params._replacementRate) {
                    // Replace weight entirely
                    attributes.weight = static_cast<float>(weightDist(gen));
                    shouldMutate = true;
                }
                break;
                
            case WeightMutationParams::MutationType::MIXED:
                // First check for replacement (lower probability)
                if (uniformDist(gen) < params._replacementRate) {
                    attributes.weight = static_cast<float>(weightDist(gen));
                    shouldMutate = true;
                } else if (uniformDist(gen) < params._perturbationRate) {
                    // Otherwise, try perturbation
                    attributes.weight += static_cast<float>(gaussianDist(gen));
                    shouldMutate = true;
                }
                break;
        }
        
        // Apply the mutation if one occurred
        if (shouldMutate) {
            // Update the connection gene with modified attributes
            // Note: This requires friend access or a setter method
            auto& mutableAttributes = connection.get_attributes();
            mutableAttributes = attributes;
            
            // Add connection history ID to deltas
            connectionDeltas.push_back(connection.get_historyID());
        }
    }
    
    // // Debug: Log final mutation results
    // static auto logger = LOGGER("operator.WeightMutation");
    // LOG_DEBUG("WeightMutation: Processed {} connections, {} mutations made (deltas: {})", 
    //          connectionGenes.size(), connectionDeltas.size(), connectionDeltas.size());
    
    return mutatedGenome;
}