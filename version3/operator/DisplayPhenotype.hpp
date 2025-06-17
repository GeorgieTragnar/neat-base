#pragma once

#include <sstream>

#include "version3/data/Genome.hpp"

namespace Operator {

void displayPhenotype(const Genome& genome, std::stringstream& output);

}