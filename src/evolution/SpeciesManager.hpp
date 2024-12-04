// species.hpp
#pragma once

#include <cstdint>
#include "core/Species.hpp"

namespace neat {
namespace core {
class Population;
class Genome;
}
namespace evolution {

class SpeciesManager {
public:
    struct Config {
        double compatibilityThreshold;
        double targetSpeciesCount;
        int32_t stagnationThreshold;
        double minFitnessThreshold;
        double activationDiffWeight;

        Config() = delete;
        Config(double compatibilityThreshold,
                double targetSpeciesCount,
                int32_t stagnationThreshold,
                double minFitnessThreshold,
                double activationDiffWeight)
            : compatibilityThreshold(compatibilityThreshold)
            , targetSpeciesCount(targetSpeciesCount)
            , stagnationThreshold(stagnationThreshold)
            , minFitnessThreshold(minFitnessThreshold)
            , activationDiffWeight(activationDiffWeight) {}
    };

    explicit SpeciesManager(const Config& config);
    
    void speciate(core::Population& population);
    void updateSpecies(core::Population& population);
    void adjustFitness(core::Population& population);
    void removeStagnantSpecies(core::Population& population);
    
    double getCompatibilityDistance(const core::Genome& genome1, const core::Genome& genome2) const;
    const std::vector<core::Species>& getSpecies() const;

private:
    void adjustCompatibilityThreshold();
    void assignGenomeToSpecies(core::Genome& genome, core::Population& population);

    Config config;
    std::vector<core::Species> species;
    double currentCompatibilityThreshold;
};

}
}
