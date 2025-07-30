#pragma once

#include <memory>

#include "version3/data/Genome.hpp"

namespace Operator {

class RepairOperatorParams;

Genome repair(const Genome& genome, const RepairOperatorParams& params);

class RepairOperatorParams {
public:
    RepairOperatorParams() = delete;
    RepairOperatorParams(uint32_t maxRepairAttempts);

protected:
    friend Genome repair(const Genome& genome, const RepairOperatorParams& params);
    
    const uint32_t _maxRepairAttempts;  // Maximum number of repair attempts before marking for elimination
};

}