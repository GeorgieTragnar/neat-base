
#include "NodeGene.hpp"

NodeGene::NodeGene(const uint32_t historyID, const NodeType type, const NodeGeneAttributes attributes)
	: _historyID(historyID)
	, _type(type)
	, _attributes(attributes)
{

}

const uint32_t& NodeGene::get_historyID() const
{
	return _historyID;
}

const NodeType& NodeGene::get_type() const
{
	return _type;
}

const NodeGeneAttributes& NodeGene::get_attributes() const
{
	return _attributes;
}

NodeGeneAttributes& NodeGene::get_attributes()
{
	return _attributes;
}
