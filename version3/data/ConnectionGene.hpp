// ConnectionGene.h
#pragma once
#include <cstdint>
#include <vector>
#include <memory>

#include "NodeGene.hpp"

#include "operator_forward_declarations.inc"

struct ConnectionGeneAttributes {
	float weight;
	bool enabled;
};

class ConnectionGene {
public:
	ConnectionGene() = delete;
	
	// Simple data container constructor
	ConnectionGene(const uint32_t historyID, size_t sourceNodeIndex, size_t targetNodeIndex, 
				const ConnectionGeneAttributes attributes);
	
	// Default copy/move operations (simple data container)
	ConnectionGene(const ConnectionGene& other) = default;
	ConnectionGene& operator=(const ConnectionGene& other) = default;
	ConnectionGene(ConnectionGene&& other) = default;
	ConnectionGene& operator=(ConnectionGene&& other) = default;
	
	const uint32_t& get_historyID() const;
	
	// Index-based accessors
	size_t get_sourceNodeIndex() const;
	size_t get_targetNodeIndex() const;
	
	const ConnectionGeneAttributes& get_attributes() const;

	bool operator==(const ConnectionGene& other) const;
	bool operator!=(const ConnectionGene& other) const;

protected:
	friend class Genome;
#include "operator_friend_declarations.inc"

	ConnectionGeneAttributes& get_attributes();

private:
	uint32_t _historyID;
	size_t _sourceNodeIndex;
	size_t _targetNodeIndex;
	ConnectionGeneAttributes _attributes;
};
