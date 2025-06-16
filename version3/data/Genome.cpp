
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
		
		const NodeGene& sourceNode = _nodeGenes[sourceIt->second];
		const NodeGene& targetNode = _nodeGenes[targetIt->second];
		
		assert(sourceNode.get_type() != NodeType::OUTPUT && "Output nodes cannot be connection sources");
		assert(targetNode.get_type() != NodeType::INPUT && "Input nodes cannot be connection targets");
		
		_connectionGenes.emplace_back(
			params._connectionHistoryIDs[i],
			sourceNode,
			targetNode,
			params._connectionAttributes[i]
		);
	}
}

Genome::Genome(const RawGenomeParams& params)
{
	assert(!params._nodeGenes.empty() && "Raw node gene data cannot be empty");
	for (const void* nodeData : params._nodeGenes) {
		assert(nodeData != nullptr && "Raw node data pointer cannot be null");
	}

	for (const void* connData : params._rawConnectionGeneData) {
		assert(connData != nullptr && "Raw connection data pointer cannot be null");
	}

	// Strategic pre-allocation to prevent reallocations during evolution
	size_t nodeCapacity = calculateOptimalCapacity(params._nodeGenes.size(), 0);
	size_t connCapacity = calculateOptimalCapacity(params._rawConnectionGeneData.size(), 0);
	
	_nodeGenes.reserve(nodeCapacity);
	_connectionGenes.reserve(connCapacity);
	
	for (const NodeGene* nodeGene : params._nodeGenes) {
		_nodeGenes.emplace_back(*nodeGene);
	}

	std::unordered_map<uint32_t, const NodeGene*> nodeMap;
	for (const NodeGene& node : _nodeGenes) {
		const uint32_t historyID = node.get_historyID();

		assert(nodeMap.find(historyID) == nodeMap.end() && 
			"Duplicate node history ID found in raw data");

		nodeMap[historyID] = &node;
	}

	for (const void* connData : params._rawConnectionGeneData) {
		const uint32_t* sourceNodeIDPtr = ConnectionGene::getSourceNodeHistoryID(connData);
		const uint32_t* targetNodeIDPtr = ConnectionGene::getTargetNodeHistoryID(connData);

		auto sourceIt = nodeMap.find(*sourceNodeIDPtr);
		auto targetIt = nodeMap.find(*targetNodeIDPtr);

		assert(sourceIt != nodeMap.end() && 
			"Connection references nonexistent source node");
		assert(targetIt != nodeMap.end() && 
			"Connection references nonexistent target node");

		const NodeGene& sourceNode = *sourceIt->second;
		const NodeGene& targetNode = *targetIt->second;

		assert(sourceNode.get_historyID() == *sourceNodeIDPtr && 
			"Source node history ID mismatch");
		assert(targetNode.get_historyID() == *targetNodeIDPtr && 
			"Target node history ID mismatch");

		assert(sourceNode.get_type() != NodeType::OUTPUT && 
			"Output nodes cannot be connection sources");
		assert(targetNode.get_type() != NodeType::INPUT && 
			"Input nodes cannot be connection targets");

		_connectionGenes.emplace_back(connData, sourceNode, targetNode);
	}
}

Genome::Genome(const Genome& other)
{
	assert(other._connectionGeneDeltas.empty() && "Source genome must have empty connection deltas for copy construction");
	assert(other._nodeGeneDeltas.empty() && "Source genome must have empty node deltas for copy construction");
	
	// Strategic pre-allocation to prevent reallocations during evolution
	size_t nodeCapacity = calculateOptimalCapacity(other._nodeGenes.size(), 0);
	size_t connCapacity = calculateOptimalCapacity(other._connectionGenes.size(), 0);
	
	_nodeGenes.reserve(nodeCapacity);
	_connectionGenes.reserve(connCapacity);
	_phenotype = other._phenotype;

	for (const NodeGene& node : other._nodeGenes) {
		_nodeGenes.emplace_back(node);
	}

	std::unordered_map<const NodeGene*, const NodeGene*> nodeMap;
	for (size_t i = 0; i < other._nodeGenes.size(); ++i) {
		LOG_DEBUG("Mapping node [{}]: old ptr {} -> new ptr {}", i, 
			static_cast<const void*>(&other._nodeGenes[i]), static_cast<const void*>(&_nodeGenes[i]));
		nodeMap[&other._nodeGenes[i]] = &_nodeGenes[i];
	}

	LOG_DEBUG("Processing {} connection genes", other._connectionGenes.size());
	for (size_t connIdx = 0; connIdx < other._connectionGenes.size(); ++connIdx) {
		const ConnectionGene& conn = other._connectionGenes[connIdx];
		const void* rawConnData = conn.get_rawData();

		LOG_DEBUG("Connection [{}]: Processing connection with raw data ptr {}", connIdx, rawConnData);
		
		const NodeGene* originalSource = &conn.get_sourceNodeGene();
		const NodeGene* originalTarget = &conn.get_targetNodeGene();
		
		LOG_DEBUG("Connection [{}]: Source ptr {}, Target ptr {}", connIdx, 
			static_cast<const void*>(originalSource), static_cast<const void*>(originalTarget));

		auto sourceIt = nodeMap.find(originalSource);
		auto targetIt = nodeMap.find(originalTarget);

		if (sourceIt == nodeMap.end() || targetIt == nodeMap.end()) {
			LOG_ERROR("Genome copy constructor node mapping failure");
			if (originalSource) {
				LOG_ERROR("Source node ID: {}, Source ptr: {}, Found: {}", 
					originalSource->get_historyID(), static_cast<const void*>(originalSource), sourceIt != nodeMap.end());
			} else {
				LOG_ERROR("Source node is NULL pointer");
			}
			if (originalTarget) {
				LOG_ERROR("Target node ID: {}, Target ptr: {}, Found: {}", 
					originalTarget->get_historyID(), static_cast<const void*>(originalTarget), targetIt != nodeMap.end());
			} else {
				LOG_ERROR("Target node is NULL pointer");
			}
			LOG_ERROR("NodeMap size: {}, Other nodes: {}, This nodes: {}", nodeMap.size(), other._nodeGenes.size(), _nodeGenes.size());
			LOG_ERROR("Available nodes in source genome:");
			for (size_t i = 0; i < other._nodeGenes.size(); ++i) {
				LOG_ERROR("  [{}] Node ID: {}, ptr: {}", i, other._nodeGenes[i].get_historyID(), static_cast<const void*>(&other._nodeGenes[i]));
			}
		}

		assert(sourceIt != nodeMap.end() && "Failed to find source node in node map");
		assert(targetIt != nodeMap.end() && "Failed to find target node in node map");

		const NodeGene& newSource = *sourceIt->second;
		const NodeGene& newTarget = *targetIt->second;

		assert(newSource.get_historyID() == originalSource->get_historyID() && 
			"Source node history ID mismatch in copy");
		assert(newTarget.get_historyID() == originalTarget->get_historyID() && 
			"Target node history ID mismatch in copy");

		_connectionGenes.emplace_back(rawConnData, *sourceIt->second, *targetIt->second);
	}

	assert(_nodeGenes.size() == other._nodeGenes.size() && 
		"Node count mismatch after copy");
	assert(_connectionGenes.size() == other._connectionGenes.size() && 
		"Connection count mismatch after copy");
}

Genome& Genome::operator=(const Genome& other)
{
	if (this == &other) {
		return *this;
	}

	assert(other._connectionGeneDeltas.empty() && "Source genome must have empty connection deltas for copy assignment");
	assert(other._nodeGeneDeltas.empty() && "Source genome must have empty node deltas for copy assignment");

	_nodeGenes.clear();
	_connectionGenes.clear();
	_connectionGeneDeltas.clear();
	_nodeGeneDeltas.clear();
	_phenotype = other._phenotype;

	// Strategic pre-allocation to prevent reallocations during evolution
	size_t nodeCapacity = calculateOptimalCapacity(other._nodeGenes.size(), 0);
	size_t connCapacity = calculateOptimalCapacity(other._connectionGenes.size(), 0);
	
	_nodeGenes.reserve(nodeCapacity);
	_connectionGenes.reserve(connCapacity);

	for (const NodeGene& node : other._nodeGenes) {
		_nodeGenes.emplace_back(node);
	}

	std::unordered_map<uint32_t, const NodeGene*> nodeHistoryMap;
	for (const NodeGene& node : _nodeGenes) {
		nodeHistoryMap[node.get_historyID()] = &node;
	}

	for (const ConnectionGene& conn : other._connectionGenes) {
		const void* rawConnData = conn.get_rawData();
		const uint32_t* sourceNodeIDPtr = ConnectionGene::getSourceNodeHistoryID(rawConnData);
		const uint32_t* targetNodeIDPtr = ConnectionGene::getTargetNodeHistoryID(rawConnData);

		auto sourceIt = nodeHistoryMap.find(*sourceNodeIDPtr);
		auto targetIt = nodeHistoryMap.find(*targetNodeIDPtr);

		assert(sourceIt != nodeHistoryMap.end() && "Failed to find source node by history ID");
		assert(targetIt != nodeHistoryMap.end() && "Failed to find target node by history ID");

		_connectionGenes.emplace_back(rawConnData, *sourceIt->second, *targetIt->second);
	}

	return *this;
}

Genome::Genome(Genome&& other) noexcept
	: _nodeGenes(std::move(other._nodeGenes))
	, _connectionGenes(std::move(other._connectionGenes))
	, _phenotype(std::move(other._phenotype))
	, _connectionGeneDeltas(std::move(other._connectionGeneDeltas))
	, _nodeGeneDeltas(std::move(other._nodeGeneDeltas))
{
}

Genome& Genome::operator=(Genome&& other) noexcept
{
	if (this == &other) {
		return *this;
	}
	
	_nodeGenes = std::move(other._nodeGenes);
	_connectionGenes = std::move(other._connectionGenes);
	_phenotype = std::move(other._phenotype);
	_connectionGeneDeltas = std::move(other._connectionGeneDeltas);
	_nodeGeneDeltas = std::move(other._nodeGeneDeltas);
	
	return *this;
}

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
