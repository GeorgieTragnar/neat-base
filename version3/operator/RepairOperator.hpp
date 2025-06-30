#pragma once

#include <memory>

#include "version3/data/Genome.hpp"
#include "version3/population/PopulationData.hpp"
#include "version3/population/GlobalIndexRegistry.hpp"

namespace Operator {

class RepairOperatorParams;

Genome repair(const Genome& genome, Population::DynamicGenomeData& genomeData, const RepairOperatorParams& params, Population::GlobalIndexRegistry& registry, uint32_t globalIndex);

class RepairOperatorParams {
public:
    RepairOperatorParams() = delete;
    RepairOperatorParams(uint32_t maxRepairAttempts);

protected:
    friend Genome repair(const Genome& genome, Population::DynamicGenomeData& genomeData, const RepairOperatorParams& params, Population::GlobalIndexRegistry& registry, uint32_t globalIndex);
    
    const uint32_t _maxRepairAttempts;  // Maximum number of repair attempts before marking for elimination
};

}