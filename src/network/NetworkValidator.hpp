
#pragma once

namespace neat {
namespace core {
class Genome;
}
namespace network {

class NetworkValidator {
public:
	// Main validation method
	static void validateNetwork(const core::Genome& genome);

private:
	// Validate all nodes exist and have valid IDs
	static void validateNodeStructure(const core::Genome& genome);

	// Validate all genes reference valid nodes and have valid structure
	static void validateGeneStructure(const core::Genome& genome);

	// Validate network topology (no cycles in feedforward networks)
	static void validateTopology(const core::Genome& genome);

	// Validate all output nodes are reachable from inputs
	static void validateConnectivity(const core::Genome& genome);
};

}
}
