// Genome.h
#pragma once
#include <vector>
#include <memory>

#include "NodeGene.hpp"
#include "ConnectionGene.hpp"
#include "Phenotype.hpp"

// Auto-generated operator forward declarations
#include "data_forward_declarations.inc"
#include "operator_forward_declarations.inc"

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
	std::vector<const NodeGene*>		_nodeGenes;
	std::vector<const ConnectionGene*>	_connectionGenes;
};

class Genome {
public:
	Genome() = delete;
	Genome(const GenomeParams& params);
	Genome(const RawGenomeParams& params);

	// Copy/move operations
	Genome(const Genome& other) = default;
	Genome& operator=(const Genome& other) = default;
	Genome(Genome&& other) noexcept = default;
	Genome& operator=(Genome&& other) noexcept = default;

	// TODO: Implement when project-wide serialization is addressed
	static Genome Deserialize(const std::vector<uint8_t>& serializedData) = delete;
	static const std::vector<uint8_t>& Serialize(const Genome& other) = delete;

	// Phenotype class moved to independent Phenotype.hpp file
	// Manual friend relationship maintained through friend declaration in Phenotype class

	const std::vector<NodeGene>& get_nodeGenes() const;
	const std::vector<ConnectionGene>& get_connectionGenes() const;
	const Phenotype& get_phenotype() const;

protected:
#include "operator_friend_declarations.inc"

	std::vector<NodeGene>& get_nodeGenes();
	std::vector<ConnectionGene>& get_connectionGenes();
	Phenotype& get_phenotype();
	
	std::vector<uint32_t>& get_connectionGeneDeltas();
	std::vector<uint32_t>& get_nodeGeneDeltas();
	
	// Capacity management for preventing vector reallocations
	void ensureCapacity(size_t additionalNodes = 0, size_t additionalConnections = 0);

private:
	// Capacity management constants
	static constexpr size_t INITIAL_NODE_CAPACITY = 64;
	static constexpr size_t INITIAL_CONN_CAPACITY = 256;
	static constexpr size_t GROWTH_MULTIPLIER = 2;
	static constexpr size_t SAFETY_MARGIN = 8;
	
	static size_t nextPowerOfTwo(size_t n);
	static size_t calculateOptimalCapacity(size_t currentSize, size_t requestedAdditional);

private:
	std::vector<NodeGene>			_nodeGenes;
	std::vector<ConnectionGene>		_connectionGenes;
	Phenotype						_phenotype;
	std::vector<uint32_t>			_connectionGeneDeltas;
	std::vector<uint32_t>			_nodeGeneDeltas;
};
