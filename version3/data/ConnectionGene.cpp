
#include "ConnectionGene.hpp"

ConnectionGene::ConnectionGene(const uint32_t historyID, const NodeGene& sourceNodeGene,
	const NodeGene& targetNodeGene, const ConnectionGeneAttributes attributes)
	: _historyID(historyID)
	, _sourceNodeGene(sourceNodeGene)
	. _targetNodeGene(targetNodeGene)
{
	
}

const uint32_t& ConnectionGene::get_historyID() const
{
	return _historyID;
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
