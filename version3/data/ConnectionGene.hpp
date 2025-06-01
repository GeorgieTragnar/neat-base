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
	ConnectionGene(const uint32_t historyID, const NodeGene& sourceNodeGene, const NodeGene& targetNodeGene, 
				const ConnectionGeneAttributes attributes);
	// void* data constructor only for runtime evolution optimization purposes
	// never use on crossplatform data
	// INTERNAL USE ONLY
	ConnectionGene(const void* data, const NodeGene& sourceNodeGene, const NodeGene& targetNodeGene);
	
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
#include "operator_friend_declarations.inc"

	static const uint32_t* getSourceNodeHistoryID(const void* data);
	static const uint32_t* getTargetNodeHistoryID(const void* data);

	const void* get_rawData() const;
	ConnectionGeneAttributes& get_attributes();

private:
	struct RawData {
		const uint32_t _historyID;
		const uint32_t _sourceNodeGeneHistoryID;
		const uint32_t _targetNodeGeneHistoryID;
		ConnectionGeneAttributes _attributes;
	} _rawData;
	const NodeGene& _sourceNodeGene;
	const NodeGene& _targetNodeGene;
};
