
#include "HistoryTracker.hpp"


uint32_t HistoryTracker::register_newNode() {
	uint32_t id = _nextNodeID++;
	return id;
}

uint32_t HistoryTracker::register_newConnection() {
	return _nextConnectionID++;
}

uint32_t HistoryTracker::get_input(size_t index) {
	if (index >= _inputIDs.size()) {
		_inputIDs.resize(index + 1, 0);
	}
	uint32_t& id = _inputIDs[index];
	return id ? id : (id = register_newNode());
}

uint32_t HistoryTracker::get_output(size_t index) {
	if (index >= _outputIDs.size()) {
		_outputIDs.resize(index + 1, 0);
	}
	uint32_t& id = _outputIDs[index];
	return id ? id : (id = register_newNode());
}

uint32_t HistoryTracker::get_bias() {
	return _biasID ? _biasID : (_biasID = register_newNode());
}

uint32_t HistoryTracker::get_connection(uint32_t sourceNodeID, uint32_t targetNodeID) {
	auto key = std::make_pair(sourceNodeID, targetNodeID);
	auto it = _connectionIDs.find(key);
	if (it != _connectionIDs.end()) return it->second;
	
	uint32_t id = register_newConnection();
	_connectionIDs[key] = id;
	return id;
}

uint32_t HistoryTracker::get_splitNode(uint32_t connectionID) {
	auto it = _splitNodeIDs.find(connectionID);
	if (it != _splitNodeIDs.end()) return it->second;

	uint32_t id = register_newNode();
	_splitNodeIDs[connectionID] = id;
	return id;
}

std::vector<uint32_t> HistoryTracker::get_allSplitNodes(uint32_t connectionID) {
	std::vector<uint32_t> result;

	auto primaryIt = _splitNodeIDs.find(connectionID);
	if (primaryIt != _splitNodeIDs.end()) {
		result.push_back(primaryIt->second);
	}
	
	auto branchIt = _splitNodeBranches.find(connectionID);
	if (branchIt != _splitNodeBranches.end()) {
		result.insert(result.end(), branchIt->second.begin(), branchIt->second.end());
	}
	
	return result;
}

uint32_t HistoryTracker::create_splitBranch(uint32_t connectionID) {
	uint32_t id = register_newNode();
	_splitNodeBranches[connectionID].push_back(id);
	return id;
}
