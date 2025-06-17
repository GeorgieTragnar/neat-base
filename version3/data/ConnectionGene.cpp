
#include "ConnectionGene.hpp"

ConnectionGene::ConnectionGene(const uint32_t historyID, size_t sourceNodeIndex,
	size_t targetNodeIndex, const ConnectionGeneAttributes attributes)
	: _historyID(historyID)
	, _sourceNodeIndex(sourceNodeIndex)
	, _targetNodeIndex(targetNodeIndex)
	, _attributes(attributes)
{
	
}

const uint32_t& ConnectionGene::get_historyID() const
{
	return _historyID;
}

size_t ConnectionGene::get_sourceNodeIndex() const
{
	return _sourceNodeIndex;
}

size_t ConnectionGene::get_targetNodeIndex() const
{
	return _targetNodeIndex;
}

const ConnectionGeneAttributes& ConnectionGene::get_attributes() const
{
	return _attributes;
}

bool ConnectionGene::operator==(const ConnectionGene& other) const
{
	return other._historyID == _historyID;
}

bool ConnectionGene::operator!=(const ConnectionGene& other) const
{
	return other._historyID != _historyID;
}

ConnectionGeneAttributes& ConnectionGene::get_attributes()
{
	return _attributes;
}
