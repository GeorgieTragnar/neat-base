
#include "Genome.hpp"
#include <unordered_map>
#include <cassert>

Genome::Genome(const GenomeParams& params)
	: _phenotype(nullptr)
{
	assert(params._nodeHistoryIDs.size() == params._nodeTypes.size() && 
		params._nodeHistoryIDs.size() == params._nodeAttributes.size() && 
		"Node parameter arrays must be of equal size");
	
	assert(params._connectionHistoryIDs.size() == params._sourceNodeHistoryIDs.size() && 
		params._connectionHistoryIDs.size() == params._targetNodeHistoryIDs.size() && 
		params._connectionHistoryIDs.size() == params._connectionAttributes.size() && 
		"Connection parameter arrays must be of equal size");

	_nodeGenes.reserve(params._nodeHistoryIDs.size());
	_connectionGenes.reserve(params._connectionHistoryIDs.size());

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
	: _phenotype(nullptr)
{
	assert(!params._rawNodeGeneData.empty() && "Raw node gene data cannot be empty");
	for (const void* nodeData : params._rawNodeGeneData) {
		assert(nodeData != nullptr && "Raw node data pointer cannot be null");
	}

	for (const void* connData : params._rawConnectionGeneData) {
		assert(connData != nullptr && "Raw connection data pointer cannot be null");
	}

	_nodeGenes.reserve(params._rawNodeGeneData.size());
	_connectionGenes.reserve(params._rawConnectionGeneData.size());
	
	for (const void* nodeData : params._rawNodeGeneData) {
		_nodeGenes.emplace_back(nodeData);
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
	: _phenotype(nullptr)
{
	_nodeGenes.reserve(other._nodeGenes.size());
	_connectionGenes.reserve(other._connectionGenes.size());

	for (const NodeGene& node : other._nodeGenes) {
		_nodeGenes.emplace_back(node.get_rawData());
	}

	std::unordered_map<const NodeGene*, const NodeGene*> nodeMap;
	for (size_t i = 0; i < other._nodeGenes.size(); ++i) {
		nodeMap[&other._nodeGenes[i]] = &_nodeGenes[i];
	}

	for (const ConnectionGene& conn : other._connectionGenes) {
		const void* rawConnData = conn.get_rawData();

		const NodeGene* originalSource = &conn.get_sourceNodeGene();
		const NodeGene* originalTarget = &conn.get_targetNodeGene();

		auto sourceIt = nodeMap.find(originalSource);
		auto targetIt = nodeMap.find(originalTarget);

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

	_nodeGenes.clear();
	_connectionGenes.clear();
	_phenotype = nullptr;

	_nodeGenes.reserve(other._nodeGenes.size());
	_connectionGenes.reserve(other._connectionGenes.size());

	for (const NodeGene& node : other._nodeGenes) {
		_nodeGenes.emplace_back(node.get_rawData());
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

const std::vector<NodeGene>& Genome::get_nodeGenes() const
{
	return _nodeGenes;
}

const std::vector<ConnectionGene>& Genome::get_connectionGenes() const
{
	return _connectionGenes;
}

std::shared_ptr<const Phenotype> Genome::get_phenotype() const
{
	return std::const_pointer_cast<const Phenotype>(_phenotype);
}

std::vector<NodeGene>& Genome::get_nodeGenes()
{
	if (_phenotype)
		_phenotype->_dirty = true;
	return _nodeGenes;
}

std::vector<ConnectionGene>& Genome::get_connectionGenes()
{
	if (_phenotype)
		_phenotype->_dirty = true;
	return _connectionGenes;
}
