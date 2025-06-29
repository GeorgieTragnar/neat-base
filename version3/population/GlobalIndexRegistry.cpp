#include "GlobalIndexRegistry.hpp"
#include <cassert>

namespace Population {

GlobalIndexRegistry::GlobalIndexRegistry(uint32_t maxIndices) 
    : _states(maxIndices, GenomeState::Active) {
}

GenomeState GlobalIndexRegistry::getState(uint32_t globalIndex) const {
    assert(globalIndex < _states.size() && "Global index out of range");
    return _states[globalIndex];
}

void GlobalIndexRegistry::markForElimination(uint32_t globalIndex) {
    assert(globalIndex < _states.size() && "Global index out of range");
    assert(_states[globalIndex] == GenomeState::Active && "Can only mark Active genomes for elimination");
    _states[globalIndex] = GenomeState::HotElimination;
}

void GlobalIndexRegistry::transitionToCold(uint32_t globalIndex) {
    assert(globalIndex < _states.size() && "Global index out of range");
    assert(_states[globalIndex] == GenomeState::HotElimination && "Can only transition from HotElimination to ColdElimination");
    _states[globalIndex] = GenomeState::ColdElimination;
}

void GlobalIndexRegistry::markReadyForReplacement(uint32_t globalIndex) {
    assert(globalIndex < _states.size() && "Global index out of range");
    assert(_states[globalIndex] == GenomeState::ColdElimination && "Can only mark ColdElimination genomes as ready for replacement");
    _states[globalIndex] = GenomeState::ReadyForReplacement;
}

void GlobalIndexRegistry::resetToActive(uint32_t globalIndex) {
    assert(globalIndex < _states.size() && "Global index out of range");
    _states[globalIndex] = GenomeState::Active;
}

uint32_t GlobalIndexRegistry::getFreeIndex() {
    // Find first ReadyForReplacement index
    for (uint32_t i = 0; i < _states.size(); ++i) {
        if (_states[i] == GenomeState::ReadyForReplacement) {
            _states[i] = GenomeState::Active;  // Reserve it immediately
            return i;
        }
    }
    return INVALID_INDEX;
}

void GlobalIndexRegistry::incrementMaxIndex() {
    _states.emplace_back(GenomeState::Active);
}

} // namespace Population