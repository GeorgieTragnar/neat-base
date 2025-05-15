// Genome.h
#pragma once
#include <vector>

#include "NodeGene.hpp"
#include "ConnectionGene.hpp"

struct GenomeParams {
	std::vector<uint32_t> 					_nodeHistoryIDs;
	std::vector<NodeType> 					_nodeTypes;
	std::vector<NodeGeneAttributes> 		_nodeAttributes;

	std::vector<uint32_t> 					_connectionHistoryIDs;
	std::vector<uint32_t> 					_sourceNodeHistoryIDs;
	std::vector<uint32_t> 					_targetNodeHistoryIDs;
	std::vector<ConnectionGeneAttributes> 	_connectionAttributes;
};

struct RawGenomeParams {
	std::vector<const void*>	_rawNodeGeneData;
	std::vector<const void*>	_rawConnectionGeneData;
};

class Genome {
public:
	Genome() = delete;
	Genome(const GenomeParams& params);
	Genome(const RawGenomeParams& params);

	Genome(const Genome& other);
	Genome& operator=(const Genome& other);
	Genome(Genome&& other) = default;
	Genome& operator=(Genome&& other) = default;

	// TODO: Implement when project-wide serialization is addressed
	static Genome Deserialize(const std::vector<uint8_t>& serializedData) = delete;
	static const std::vector<uint8_t>& Serialize(const Genome& other) = delete;


private:
	std::vector<NodeGene> _nodeGenes;
	std::vector<ConnectionGene> _connectionGenes;
};
