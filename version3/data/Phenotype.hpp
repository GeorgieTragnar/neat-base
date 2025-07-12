#pragma once

#include <vector>
#include <cstddef>

// Forward declarations
class Genome;

// Need full definitions for struct fields
#include "ConnectionGene.hpp"  // For ConnectionGeneAttributes
#include "NodeGene.hpp"        // For NodeGeneAttributes

class Phenotype {
public:	
	struct Connection {
		size_t						_sourceNodeIndex;
		size_t						_targetNodeIndex;
		ConnectionGeneAttributes	_connectionGeneAttribute;
	};
	std::vector<NodeGeneAttributes>	_nodeGeneAttributes;
	std::vector<Connection> 		_orderedConnections;

	std::vector<size_t>				_inputIndices;
	std::vector<size_t>				_outputIndices;
	size_t							_biasIndex;

private:
	// Manual friend declaration to maintain the "view" relationship
	// Genome can construct and access Phenotype as its compiled representation
	friend class Genome;
};