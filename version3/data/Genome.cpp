
#include "Genome.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cassert>

#include "../logger/Logger.hpp"
static auto logger = LOGGER("core::Genome");

Genome::Genome(const GenomeParams& params)
{
	assert(params._nodeHistoryIDs.size() == params._nodeTypes.size() && 
		params._nodeHistoryIDs.size() == params._nodeAttributes.size() && 
		"Node parameter arrays must be of equal size");
	
	assert(params._connectionHistoryIDs.size() == params._sourceNodeHistoryIDs.size() && 
		params._connectionHistoryIDs.size() == params._targetNodeHistoryIDs.size() && 
		params._connectionHistoryIDs.size() == params._connectionAttributes.size() && 
		"Connection parameter arrays must be of equal size");

	// Strategic pre-allocation to prevent reallocations during evolution
	size_t nodeCapacity = calculateOptimalCapacity(params._nodeHistoryIDs.size(), 0);
	size_t connCapacity = calculateOptimalCapacity(params._connectionHistoryIDs.size(), 0);
	
	_nodeGenes.reserve(nodeCapacity);
	_connectionGenes.reserve(connCapacity);

	for (size_t i = 0; i < params._nodeHistoryIDs.size(); ++i) {
		_nodeGenes.emplace_back(
			params._nodeHistoryIDs[i],
			params._nodeTypes[i],
			params._nodeAttributes[i]
		);
	}

	std::unordered_map<uint32_t, size_t> nodeIndexMap;
	for (size_t i = 0; i < _nodeGenes.size(); ++i) {
		nodeIndexMap[_nodeGenes[i].get_historyID()] = i;
	}

	for (size_t i = 0; i < params._connectionHistoryIDs.size(); ++i) {
		uint32_t sourceNodeID = params._sourceNodeHistoryIDs[i];
		uint32_t targetNodeID = params._targetNodeHistoryIDs[i];
		
		auto sourceIt = nodeIndexMap.find(sourceNodeID);
		auto targetIt = nodeIndexMap.find(targetNodeID);
		
		assert(sourceIt != nodeIndexMap.end() && "Connection references nonexistent source node");
		assert(targetIt != nodeIndexMap.end() && "Connection references nonexistent target node");
		
		size_t sourceIndex = sourceIt->second;
		size_t targetIndex = targetIt->second;
		
		const NodeGene& sourceNode = _nodeGenes[sourceIndex];
		const NodeGene& targetNode = _nodeGenes[targetIndex];
		
		assert(sourceNode.get_type() != NodeType::OUTPUT && "Output nodes cannot be connection sources");
		assert(targetNode.get_type() != NodeType::INPUT && "Input nodes cannot be connection targets");
		
		_connectionGenes.emplace_back(
			params._connectionHistoryIDs[i],
			sourceIndex,
			targetIndex,
			params._connectionAttributes[i]
		);
	}
}

Genome::Genome(const RawGenomeParams& params)
{
	assert(!params._nodeGenes.empty() && "Raw node gene data cannot be empty");
	for (const NodeGene* nodePtr : params._nodeGenes) {
		assert(nodePtr != nullptr && "Raw node data pointer cannot be null");
	}

	for (const ConnectionGene* connPtr : params._connectionGenes) {
		assert(connPtr != nullptr && "Raw connection data pointer cannot be null");
	}

	// Strategic pre-allocation to prevent reallocations during evolution
	size_t nodeCapacity = calculateOptimalCapacity(params._nodeGenes.size(), 0);
	size_t connCapacity = calculateOptimalCapacity(params._connectionGenes.size(), 0);
	
	_nodeGenes.reserve(nodeCapacity);
	_connectionGenes.reserve(connCapacity);
	
	// Simple copying - no complex mapping needed with index-based design
	for (const NodeGene* nodePtr : params._nodeGenes) {
		_nodeGenes.emplace_back(*nodePtr);
	}

	// Check for duplicate node history IDs
	std::unordered_set<uint32_t> seenHistoryIDs;
	for (const NodeGene& node : _nodeGenes) {
		assert(seenHistoryIDs.find(node.get_historyID()) == seenHistoryIDs.end() && 
			"Duplicate node history ID found in raw data");
		seenHistoryIDs.insert(node.get_historyID());
	}

	for (const ConnectionGene* connPtr : params._connectionGenes) {
		_connectionGenes.emplace_back(*connPtr);  // Simple copy of data container
	}

	// Validation: ensure connection indices are valid for this genome
	for (const ConnectionGene& conn : _connectionGenes) {
		assert(conn.get_sourceNodeIndex() < _nodeGenes.size() && 
			"Connection source index out of bounds");
		assert(conn.get_targetNodeIndex() < _nodeGenes.size() && 
			"Connection target index out of bounds");
		
		const NodeGene& sourceNode = _nodeGenes[conn.get_sourceNodeIndex()];
		const NodeGene& targetNode = _nodeGenes[conn.get_targetNodeIndex()];
		
		assert(sourceNode.get_type() != NodeType::OUTPUT && 
			"Output nodes cannot be connection sources");
		assert(targetNode.get_type() != NodeType::INPUT && 
			"Input nodes cannot be connection targets");
	}
}

// Copy and move operations are now defaulted in header - much simpler with index-based design!

const std::vector<NodeGene>& Genome::get_nodeGenes() const
{
	return _nodeGenes;
}

const std::vector<ConnectionGene>& Genome::get_connectionGenes() const
{
	return _connectionGenes;
}

const Genome::Phenotype& Genome::get_phenotype() const
{
	return _phenotype;
}

std::vector<NodeGene>& Genome::get_nodeGenes()
{
	return _nodeGenes;
}

std::vector<ConnectionGene>& Genome::get_connectionGenes()
{
	return _connectionGenes;
}

Genome::Phenotype& Genome::get_phenotype()
{
	return _phenotype;
}

std::vector<uint32_t>& Genome::get_connectionGeneDeltas()
{
	return _connectionGeneDeltas;
}

std::vector<uint32_t>& Genome::get_nodeGeneDeltas()
{
	return _nodeGeneDeltas;
}


void Genome::ensureCapacity(size_t additionalNodes, size_t additionalConnections) {
	// Check if we need to expand node capacity
	size_t requiredNodeCapacity = _nodeGenes.size() + additionalNodes + SAFETY_MARGIN;
	if (requiredNodeCapacity > _nodeGenes.capacity()) {
		size_t newNodeCapacity = calculateOptimalCapacity(_nodeGenes.size(), additionalNodes);
		_nodeGenes.reserve(newNodeCapacity);
	}
	
	// Check if we need to expand connection capacity
	size_t requiredConnCapacity = _connectionGenes.size() + additionalConnections + SAFETY_MARGIN;
	if (requiredConnCapacity > _connectionGenes.capacity()) {
		size_t newConnCapacity = calculateOptimalCapacity(_connectionGenes.size(), additionalConnections);
		_connectionGenes.reserve(newConnCapacity);
	}
}

size_t Genome::nextPowerOfTwo(size_t n) {
	if (n == 0) return 1;
	if (n == 1) return 2;
	
	// Find the most significant bit position
	size_t msb = 0;
	size_t temp = n - 1;
	while (temp > 0) {
		temp >>= 1;
		msb++;
	}
	
	return 1ULL << msb;
}

size_t Genome::calculateOptimalCapacity(size_t currentSize, size_t requestedAdditional) {
	size_t minimumRequired = currentSize + requestedAdditional;
	
	// For small genomes, use fixed initial capacities
	size_t minNodeCapacity = std::max(minimumRequired, INITIAL_NODE_CAPACITY);
	size_t minConnCapacity = std::max(minimumRequired, INITIAL_CONN_CAPACITY);
	
	// Use appropriate minimum based on context (heuristic: connections usually outnumber nodes)
	size_t effectiveMinimum = (minimumRequired < 100) ? minConnCapacity : minNodeCapacity;
	
	// For larger genomes, use power-of-two sizing with growth multiplier
	if (minimumRequired > effectiveMinimum) {
		size_t powerOfTwoSize = nextPowerOfTwo(minimumRequired);
		// Apply growth multiplier for extra headroom
		return powerOfTwoSize * GROWTH_MULTIPLIER;
	}
	
	return effectiveMinimum;
}
