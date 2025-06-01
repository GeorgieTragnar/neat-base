// InitOperator.h
#pragma once
#include <vector>
#include <unordered_map>
#include <memory>

#include "version3/data/NodeGene.hpp"
#include "version3/data/ConnectionGene.hpp"
#include "version3/data/Genome.hpp"
#include "version3/data/HistoryTracker.hpp"

namespace Operator {

class InitializationParams;

Genome init(std::unique_ptr<HistoryTracker>& historyTracker, const InitializationParams& params);

class InitializationParams {
public:

	enum class InputConnectionStrategy {
		NONE,
		CONNECT_TO_OUTPUTS
	};

	InitializationParams() = delete;
	InitializationParams(const std::vector<NodeGeneAttributes>& inputAttributes, 
		const std::vector<NodeGeneAttributes>& outputAttributes,
		const std::unordered_map<size_t, ConnectionGeneAttributes>& biasAttributes,
		const InputConnectionStrategy& inputStrategy);


protected:
	friend Genome init(std::unique_ptr<HistoryTracker>& historyTracker, const InitializationParams& params);

	const std::vector<NodeGeneAttributes>							_inputAttributes;
	const std::vector<NodeGeneAttributes>							_outputAttributes;
	const std::unordered_map<size_t, ConnectionGeneAttributes>		_biasAttributes;

	const InputConnectionStrategy									_inputStrategy;
};

}
