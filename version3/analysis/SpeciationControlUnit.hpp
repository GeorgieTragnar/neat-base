#pragma once
#include <memory>
#include <vector>
#include "../data/Phenotype.hpp"

namespace Analysis {

class SpeciationControlUnit {
public:
    virtual ~SpeciationControlUnit() = default;
    
    virtual std::vector<std::shared_ptr<const Phenotype>> getChampions() const = 0;
    virtual std::shared_ptr<const Phenotype> getBestChampion() const = 0;
    virtual std::shared_ptr<const Phenotype> getRandomChampion() const = 0;
    virtual size_t getChampionCount() const = 0;
};

}