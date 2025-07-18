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
    assert((_states[globalIndex] == GenomeState::Active || _states[globalIndex] == GenomeState::Elite) && 
           "Can only mark Active or Elite genomes for elimination");
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

uint32_t GlobalIndexRegistry::incrementMaxIndex() {
    _states.emplace_back(GenomeState::Active);
    return static_cast<uint32_t>(_states.size() - 1);  // Return the new index
}

void GlobalIndexRegistry::markAsElite(uint32_t globalIndex) {
    assert(globalIndex < _states.size() && "Global index out of range");
    assert(_states[globalIndex] == GenomeState::Active && "Can only mark Active genomes as Elite");
    _states[globalIndex] = GenomeState::Elite;
}

void GlobalIndexRegistry::clearAllEliteStatus() {
    for (auto& state : _states) {
        if (state == GenomeState::Elite) {
            state = GenomeState::Active;
        }
    }
}

} // namespace Population