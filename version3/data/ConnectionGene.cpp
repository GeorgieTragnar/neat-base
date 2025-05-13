
#include "ConnectionGene.hpp"

ConnectionGene::ConnectionGene(const uint32_t historyID, const NodeGene& sourceNodeGene,
	const NodeGene& targetNodeGene, const ConnectionGeneAttributes attributes)
	: _rawData(
		historyID,
		sourceNodeGene.get_historyID(),
		targetNodeGene.get_historyID(),
		attributes
	)
	, _sourceNodeGene(sourceNodeGene)
	, _targetNodeGene(targetNodeGene)
{
	
}

ConnectionGene::ConnectionGene(const void* data, const NodeGene& sourceNodeGene,
	const NodeGene& targetNodeGene)
	: _rawData(*reinterpret_cast<const RawData*>(data))
	, _sourceNodeGene(sourceNodeGene)
	, _targetNodeGene(targetNodeGene)
{

}

const uint32_t* ConnectionGene::getSourceNodeHistoryID(const void* data)
{
	return &reinterpret_cast<const RawData*>(data)->_sourceNodeGeneHistoryID;
}

const uint32_t* ConnectionGene::getTargetNodeHistoryID(const void* data)
{
	return &reinterpret_cast<const RawData*>(data)->_targetNodeGeneHistoryID;
}

const uint32_t& ConnectionGene::get_historyID() const
{
	return _rawData._historyID;
}

const NodeGene& ConnectionGene::get_sourceNodeGene() const
{
	return _sourceNodeGene;
}

const NodeGene& ConnectionGene::get_targetNodeGene() const
{
	return _targetNodeGene;
}

const ConnectionGeneAttributes& ConnectionGene::get_attributes() const
{
	return _rawData._attributes;
}

bool ConnectionGene::operator==(const ConnectionGene& other) const
{
	return other._rawData._historyID == _rawData._historyID;
}

bool ConnectionGene::operator!=(const ConnectionGene& other) const
{
	return other._rawData._historyID != _rawData._historyID;
}

const void* ConnectionGene::get_rawData() const
{
	return &_rawData;
}

ConnectionGeneAttributes& ConnectionGene::get_attributes()
{
	return _rawData._attributes;
}
