// NodeGene.h
#pragma once
#include <cstdint>

enum class NodeType {
	INPUT,
	HIDDEN,
	OUTPUT,
	BIAS
};

enum class ActivationType {
	SIGMOID,
	TANH,
	RELU,
	LEAKY_RELU,
	STEP
};

class NodeGeneAttributes {
	ActivationType activationType;
};

class NodeGene {
public:
	NodeGene() = delete;
	NodeGene(const uint32_t historyID, const NodeType type, const NodeGeneAttributes attributes);
	// void* data constructor only for runtime evolution optimization purposes
	// never use on crossplatform data
	NodeGene(const void* data);
	
	NodeGene(const NodeGene& other) = default;
	NodeGene& operator=(const NodeGene&) = default;
	NodeGene(NodeGene&& other) = default;
	NodeGene& operator=(NodeGene&&) = default;
	
	const uint32_t& get_historyID() const;
	const NodeType& get_type() const;
	const NodeGeneAttributes& get_attributes() const;

protected:
	friend class Genome;

	const void* get_rawData() const;
	NodeGeneAttributes& get_attributes();

private:
	struct RawData {
		const uint32_t _historyID;
		const NodeType _type;
		NodeGeneAttributes _attributes;
	} _rawData;
};
