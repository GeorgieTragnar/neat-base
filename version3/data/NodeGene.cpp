
#include "NodeGene.hpp"

NodeGene::NodeGene(uint32_t historyID, NodeType type, NodeGeneAttributes attributes)
	: _rawData(
		historyID,
		type,
		attributes
	)
{

}

NodeGene::NodeGene(const void* data)
	: _rawData(*reinterpret_cast<const RawData*>(data))
{

}

const uint32_t& NodeGene::get_historyID() const
{
	return _rawData._historyID;
}

const NodeType& NodeGene::get_type() const
{
	return _rawData._type;
}

const NodeGeneAttributes& NodeGene::get_attributes() const
{
	return _rawData._attributes;
}

const void* NodeGene::get_rawData() const
{
	return &_rawData;
}

NodeGeneAttributes& NodeGene::get_attributes()
{
	return _rawData._attributes;
}
