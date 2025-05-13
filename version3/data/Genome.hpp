// Genome.h
#pragma once
#include <vector>

struct GenomeParams {
	const std::vector<uint32_t>	_nodeGeneHistoryIDs;
	const std::vector<uint32_t>	_connectionGeneHistoryIDs;
	const std::vector<bool>		_connectionGeneEnable;
};

class Genome {
public:
	Genome() = delete;
	Genome(const GenomeParams& params);

	Genome(const Genome& other);
	Genome& operator=(const Genome& other);
	Genome(Genome&& other) = default;
	Genome& operator=(Genome&& other) = default;


private:
	
};
