#pragma once

#include "version3/data/Genome.hpp"

namespace Operator {

// Check if two genomes are structurally identical
bool genomeEquals(const Genome& genome1, const Genome& genome2);

}