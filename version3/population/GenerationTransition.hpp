#pragma once

#include <cstdint>
#include <cassert>
#include "version3/data/GlobalIndexRegistry.hpp"

namespace Operator {

// Forward declarations
class GenerationTransitionParams;

void generationTransition(
    GlobalIndexRegistry& registry,
    size_t genomeIndex,
    const GenerationTransitionParams& params
);

// Parameters class definition
class GenerationTransitionParams {
public:
    GenerationTransitionParams() = default;
    
    GenerationTransitionParams(const GenerationTransitionParams& other) = default;
    GenerationTransitionParams& operator=(const GenerationTransitionParams& other) = default;

protected:
    // No parameters needed for basic state transitions
    
    friend void generationTransition(
        GlobalIndexRegistry& registry,
        size_t genomeIndex,
        const GenerationTransitionParams& params
    );
};

// Implementation
inline void generationTransition(
    GlobalIndexRegistry& registry,
    size_t genomeIndex,
    const GenerationTransitionParams& params
) {
    // Validate index bounds
    assert(genomeIndex < registry.getMaxIndex() && 
           "generationTransition: genomeIndex out of bounds");
    
    GenomeState currentState = registry.getState(static_cast<uint32_t>(genomeIndex));
    
    // Assert this operator is only called for non-active genomes
    assert(currentState != GenomeState::Active && 
           "generationTransition should not be called on Active genomes");
    
    switch (currentState) {
        case GenomeState::HotElimination:
            // One generation in HotElimination -> transition to ColdElimination
            registry.transitionToCold(static_cast<uint32_t>(genomeIndex));
            break;
            
        case GenomeState::ColdElimination:
            // One generation in ColdElimination -> ready for replacement
            registry.markReadyForReplacement(static_cast<uint32_t>(genomeIndex));
            break;
            
        case GenomeState::ReadyForReplacement:
            // Already ready for replacement, no further transition needed
            break;
            
        case GenomeState::Elite:
            // Elite genomes don't transition in this context
            // They remain Elite until explicitly cleared
            break;
            
        case GenomeState::Active:
            // Should not reach here due to assertion above
            assert(false && "Logic error: Active genome in generationTransition");
            break;
    }
}

} // namespace Operator