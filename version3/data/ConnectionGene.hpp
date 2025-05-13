// ConnectionGene.h
#pragma once
#include <cstdint>

#include "NodeGene.hpp"

class ConnectionGeneAttributes {
	float weight;
	bool enabled;
};

class ConnectionGene {
public:
	ConnectionGene() = delete;
	ConnectionGene(const uint32_t historyID, const NodeGene& sourceNodeGene, const NodeGene& targetNodeGene, 
				const ConnectionGeneAttributes attributes);
	
	ConnectionGene(const ConnectionGene& other) = default;
	ConnectionGene& operator=(const ConnectionGene& other) = default;
	ConnectionGene(ConnectionGene&& other) = default;
	ConnectionGene& operator=(ConnectionGene&& other) = default;
	
	const uint32_t& get_historyID() const;
	const NodeGene& get_sourceNodeGene() const;
	const NodeGene& get_targetNodeGene() const;
	
	const ConnectionGeneAttributes& get_attributes() const;

	bool operator==(const ConnectionGene& other) const;
	bool operator!=(const ConnectionGene& other) const;

protected:
	friend class Genome;

	ConnectionGeneAttributes& get_attributes();

private:
	const uint32_t _historyID;
	const NodeGene& _sourceNodeGene;
	const NodeGene& _targetNodeGene;
	ConnectionGeneAttributes _attributes;
};
