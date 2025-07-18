//HistoryTracker.h
#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>

class Genome;
class HistoryTracker;

namespace Operator {
	class CompatibilityDistanceParams;
	uint32_t compatibilityDistance(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const CompatibilityDistanceParams& params);
}

class HistoryTracker {
public:
	uint32_t get_input(size_t inputIndex);
	uint32_t get_output(size_t outputIndex);
	uint32_t get_bias();

	uint32_t get_connection(uint32_t sourceNodeID, uint32_t targetNodeID);
	uint32_t get_splitNode(uint32_t connectionID);

	std::vector<uint32_t> get_allSplitNodes(uint32_t connectionID);
	uint32_t create_splitBranch(uint32_t connectionID);

private:
	friend uint32_t Operator::compatibilityDistance(const Genome& genome, std::shared_ptr<HistoryTracker> historyTracker, const Operator::CompatibilityDistanceParams& params);

	uint32_t register_newNode();
	uint32_t register_newConnection();

	struct PairHash {
		size_t operator()(const std::pair<uint32_t, uint32_t>& p) const {
			size_t h1 = std::hash<uint32_t>{}(p.first);
			size_t h2 = std::hash<uint32_t>{}(p.second);
			return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
		}
	};

	uint32_t	_nextNodeID = 1;
	uint32_t	_nextConnectionID = 1;
	uint32_t	_nextSpeciesID = 0;

	std::vector<uint32_t>	_inputIDs;
	std::vector<uint32_t>	_outputIDs;
	uint32_t				_biasID = 0;

	std::unordered_map<std::pair<uint32_t,uint32_t>, uint32_t, PairHash>	_connectionIDs;
	std::unordered_map<uint32_t, uint32_t>									_splitNodeIDs;
	std::unordered_map<uint32_t, std::vector<uint32_t>>						_splitNodeBranches;
	
	std::vector<Genome>														_speciesRepresentatives;
};
