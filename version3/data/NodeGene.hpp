// NodeGene.h
#pragma once
#include <cstdint>
#include <memory>

#include "operator_forward_declarations.inc"

enum class NodeType {
	INPUT,
	HIDDEN,
	OUTPUT,
	BIAS
};

enum class ActivationType {
	NONE,
	SIGMOID,
	TANH,
	RELU,
	LEAKY_RELU,
	STEP
};

struct NodeGeneAttributes {
	ActivationType activationType;
};

class NodeGene {
public:
	NodeGene() = delete;
	NodeGene(const uint32_t historyID, const NodeType type, const NodeGeneAttributes attributes);
	
	NodeGene(const NodeGene& other) = default;
	NodeGene& operator=(const NodeGene& other) = default;
	NodeGene(NodeGene&& other) = default;
	NodeGene& operator=(NodeGene&& other) = default;
	
	const uint32_t& get_historyID() const;
	const NodeType& get_type() const;
	const NodeGeneAttributes& get_attributes() const;

protected:
	friend class Genome;
#include "operator_friend_declarations.inc"

	NodeGeneAttributes& get_attributes();

private:
	uint32_t _historyID;
	NodeType _type;
	NodeGeneAttributes _attributes;
};
