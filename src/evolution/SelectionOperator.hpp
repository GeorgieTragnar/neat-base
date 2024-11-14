// selection.hpp
#pragma once

#include <cstdint>
#include <vector>
#include <random>
#include <memory>

namespace neat {
namespace core {
class Genome;
class Population;
class Species;
}
namespace evolution {

class SelectionOperator {
public:
    struct Config {
        int32_t tournamentSize = 3;
        double elitismRate = 0.1;
    };

    explicit SelectionOperator(const Config& config);
    
    const core::Genome& selectParent(const core::Population& population);
    const core::Genome& tournamentSelect(const std::vector<core::Genome>& candidates);
    std::vector<std::weak_ptr<core::Genome>> selectElites(const core::Population& population);
    
    // Species selection
    int32_t selectSpecies(const std::vector<core::Species>& species);
    const core::Genome& selectFromSpecies(const core::Species& species, const core::Population& population);

private:
    Config config;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> dist;
};

}
}
