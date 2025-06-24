
#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <unordered_set>
#include <algorithm>

#include "logger/Logger.hpp"
#include <fmt/format.h>

// Core NEAT components
#include "data/Genome.hpp"
#include "data/NodeGene.hpp"
#include "data/ConnectionGene.hpp"
#include "operator/Init.hpp"
#include "operator/WeightMutation.hpp"
#include "operator/DisplayGenome.hpp"
#include "operator/DisplayPhenotype.hpp"
#include "population/DynamicDataUpdate.hpp"
#include "population/PlotElites.hpp"
#include "population/PlotCrossover.hpp"
#include "operator/CompatibilityDistance.hpp"
#include "operator/Crossover.hpp"
#include "operator/CycleDetection.hpp"
#include "operator/NodeMutation.hpp"
#include "operator/ConnectionMutation.hpp"
#include "operator/ConnectionReactivation.hpp"
#include "operator/RepairOperator.hpp"
#include "analysis/strategies/SingleSimpleFitnessStrategy.hpp"

// Evolution prototype
#include "evolution/EvolutionPrototype.hpp"

static auto logger = LOGGER("evolution.EvolutionPrototype");

// User-defined fitness result type for XOR problem
class DoubleFitnessResult {
public:
    explicit DoubleFitnessResult(double value = 0.0) : _value(value) {}
    
    bool isBetterThan(const DoubleFitnessResult& other) const {
        return _value > other._value; // Higher is better
    }
    
    bool isEqualTo(const DoubleFitnessResult& other) const {
        return std::abs(_value - other._value) < std::numeric_limits<double>::epsilon();
    }
    
    double getValue() const { return _value; }
    
    // Comparison operators for multimap ordering
    bool operator<(const DoubleFitnessResult& other) const {
        return _value < other._value;
    }

private:
    double _value;
};

// Logger configuration map
const std::map<std::string, spdlog::level::level_enum> LOG_LEVELS = {
    {"trace", spdlog::level::trace},
    {"debug", spdlog::level::debug},
    {"info", spdlog::level::info},
    {"warn", spdlog::level::warn},
    {"error", spdlog::level::err},
    {"critical", spdlog::level::critical}
};

void printUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n"
              << "Options:\n"
              << "  <logger>=<level>  Set log level for specific logger\n"
              << "Available levels: trace, debug, info, warn, error\n"
              << "Example: " << program << " main=debug core.Genome=trace\n";
}

void configureLoggers(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        std::string param = argv[i];
        
        if (param == "--help") {
            printUsage(argv[0]);
            exit(0);
        }
        
        size_t pos = param.find('=');
        if (pos != std::string::npos) {
            std::string loggerName = param.substr(0, pos);
            std::string levelStr = param.substr(pos + 1);
            auto levelIt = LOG_LEVELS.find(levelStr);
            if (levelIt != LOG_LEVELS.end()) {
                Logger::instance().set_level(loggerName, levelIt->second);
                LOG_INFO("Set log level for {} to {}", loggerName, levelStr);
            }
        }
    }
}

// XOR truth table inputs/outputs
const std::vector<std::vector<double>> XOR_INPUTS = {
    {0, 0}, {0, 1}, {1, 0}, {1, 1}
};
const std::vector<double> XOR_OUTPUTS = {0, 1, 1, 0};

// Forward declaration
double simulateNetworkOutput(const Genome::Phenotype& phenotype, const std::vector<double>& inputs);

// XOR fitness evaluation function
DoubleFitnessResult evaluateXORFitness(const Genome::Phenotype& phenotype) {
    double fitness = 0.0;
    double totalError = 0.0;
    
    // Debug: Print network outputs for each input (only for best genome)
    static bool debugPrinted = false;
    
    // Test each XOR input pattern
    for (size_t i = 0; i < XOR_INPUTS.size(); ++i) {
        const auto& inputs = XOR_INPUTS[i];
        double expectedOutput = XOR_OUTPUTS[i];
        
        double networkOutput = simulateNetworkOutput(phenotype, inputs);
        
        // Calculate squared error
        double error = networkOutput - expectedOutput;
        totalError += error * error;
    }
    
    // Convert error to fitness (lower error = higher fitness)
    double performanceFitness = 4.0 - totalError; // Max fitness is 4.0 (perfect)
    
    // Debug output for high-fitness networks (always show for debugging)
    if (performanceFitness > 2.9) {
        static int debugCount = 0;
        if (debugCount < 3) { // Show first 3 high-fitness networks
            LOG_DEBUG("=== HIGH FITNESS NETWORK #{} (fitness: {:.3f}) ===", debugCount + 1, performanceFitness);
            
            // Show network structure
            LOG_DEBUG("Network structure: {} nodes, {} connections", 
                     phenotype._nodeGeneAttributes.size(), phenotype._orderedConnections.size());
            // Manually format indices since fmt::join might not be available
            std::string inputIndicesStr = "[";
            for (size_t i = 0; i < phenotype._inputIndices.size(); ++i) {
                if (i > 0) inputIndicesStr += ", ";
                inputIndicesStr += std::to_string(phenotype._inputIndices[i]);
            }
            inputIndicesStr += "]";
            
            std::string outputIndicesStr = "[";
            for (size_t i = 0; i < phenotype._outputIndices.size(); ++i) {
                if (i > 0) outputIndicesStr += ", ";
                outputIndicesStr += std::to_string(phenotype._outputIndices[i]);
            }
            outputIndicesStr += "]";
            
            LOG_DEBUG("Input indices: {}, Output indices: {}", inputIndicesStr, outputIndicesStr);
                     
            // Show all connections
            for (size_t i = 0; i < phenotype._orderedConnections.size(); ++i) {
                const auto& conn = phenotype._orderedConnections[i];
                LOG_DEBUG("Connection[{}]: {} -> {} (weight: {:.3f}, enabled: {})",
                         i, conn._sourceNodeIndex, conn._targetNodeIndex, 
                         conn._connectionGeneAttribute.weight, conn._connectionGeneAttribute.enabled);
            }
            
            for (size_t i = 0; i < XOR_INPUTS.size(); ++i) {
                const auto& inputs = XOR_INPUTS[i];
                double expectedOutput = XOR_OUTPUTS[i];
                double networkOutput = simulateNetworkOutput(phenotype, inputs);
                LOG_DEBUG("Input: [{:.0f}, {:.0f}] → Expected: {:.0f}, Actual: {:.3f}, Error: {:.3f}", 
                         inputs[0], inputs[1], expectedOutput, networkOutput, 
                         std::abs(networkOutput - expectedOutput));
            }
            debugCount++;
        }
    }
    
    
    // Add complexity penalization
    size_t nodeCount = phenotype._nodeGeneAttributes.size();
    size_t connectionCount = 0;
    
    // Count only enabled connections
    for (const auto& conn : phenotype._orderedConnections) {
        if (conn._connectionGeneAttribute.enabled) {
            connectionCount++;
        }
    }
    
    // Very light complexity penalization (allow growth for XOR)
    double complexityPenalty = (nodeCount > 8 ? (nodeCount - 8) * 0.01 : 0.0) + 
                              (connectionCount > 10 ? (connectionCount - 10) * 0.005 : 0.0);
    
    fitness = std::max(0.0, performanceFitness - complexityPenalty); // Ensure non-negative fitness
    
    return DoubleFitnessResult(fitness);
}

// Proper feed-forward network activation with single-pass processing
double simulateNetworkOutput(const Genome::Phenotype& phenotype, const std::vector<double>& inputs) {
    const auto& connections = phenotype._orderedConnections;
    const auto& nodes = phenotype._nodeGeneAttributes;
    
    static bool detailedDebug = true; // Set to true for first network only
    
    // Initialize node values
    std::vector<double> nodeValues(nodes.size(), 0.0);
    
    // Set input values
    for (size_t i = 0; i < inputs.size() && i < phenotype._inputIndices.size(); ++i) {
        nodeValues[phenotype._inputIndices[i]] = inputs[i];
        if (detailedDebug) {
            LOG_DEBUG("Set input node[{}] = {:.3f}", phenotype._inputIndices[i], inputs[i]);
        }
    }
    
    // Find bias node: it should be the node that is not input, not output, 
    // and has outgoing connections but no incoming connections
    for (size_t i = 0; i < nodes.size(); ++i) {
        bool isInput = std::find(phenotype._inputIndices.begin(), phenotype._inputIndices.end(), i) != phenotype._inputIndices.end();
        bool isOutput = std::find(phenotype._outputIndices.begin(), phenotype._outputIndices.end(), i) != phenotype._outputIndices.end();
        
        if (!isInput && !isOutput) {
            // Check if this node has outgoing connections but no incoming connections
            bool hasOutgoingConnections = false;
            bool hasIncomingConnections = false;
            
            for (const auto& conn : connections) {
                if (conn._connectionGeneAttribute.enabled) {
                    if (conn._sourceNodeIndex == i) hasOutgoingConnections = true;
                    if (conn._targetNodeIndex == i) hasIncomingConnections = true;
                }
            }
            
            // Bias node: has outgoing connections but no incoming connections
            if (hasOutgoingConnections && !hasIncomingConnections) {
                nodeValues[i] = 1.0; // Set bias value
                if (detailedDebug) {
                    LOG_DEBUG("Found bias node[{}] = 1.0 (outgoing: {}, incoming: {})", i, hasOutgoingConnections, hasIncomingConnections);
                }
            }
        }
    }
    
    if (detailedDebug) {
        LOG_DEBUG("Initial node values after bias setup:");
        for (size_t i = 0; i < nodeValues.size(); ++i) {
            LOG_DEBUG("  node[{}] = {:.3f}", i, nodeValues[i]);
        }
    }
    
    // Create sets for faster lookup
    std::unordered_set<size_t> inputSet(phenotype._inputIndices.begin(), phenotype._inputIndices.end());
    std::unordered_set<size_t> outputSet(phenotype._outputIndices.begin(), phenotype._outputIndices.end());
    
    // Single pass: accumulate inputs for each non-input node, then apply activation
    std::vector<double> nodeInputs(nodes.size(), 0.0);
    
    // Accumulate weighted inputs for all connections
    for (const auto& conn : connections) {
        if (conn._connectionGeneAttribute.enabled) {
            size_t sourceIdx = conn._sourceNodeIndex;
            size_t targetIdx = conn._targetNodeIndex;
            
            if (sourceIdx < nodeValues.size() && targetIdx < nodeValues.size()) {
                // Skip if target is an input node (they shouldn't receive connections)
                if (inputSet.find(targetIdx) != inputSet.end()) continue;
                
                // Accumulate weighted input
                double weightedInput = nodeValues[sourceIdx] * conn._connectionGeneAttribute.weight;
                nodeInputs[targetIdx] += weightedInput;
                
                if (detailedDebug) {
                    LOG_DEBUG("Connection: node[{}]({:.3f}) -> node[{}], weight={:.3f}, contribution={:.3f}, total_input={:.3f}",
                             sourceIdx, nodeValues[sourceIdx], targetIdx, conn._connectionGeneAttribute.weight, 
                             weightedInput, nodeInputs[targetIdx]);
                }
            }
        }
    }
    
    // Apply activation functions to non-input nodes (but not bias nodes)
    for (size_t i = 0; i < nodes.size(); ++i) {
        // Skip input nodes - they keep their original values
        if (inputSet.find(i) != inputSet.end()) continue;
        
        // Check if this is a bias node (has outgoing but no incoming connections)
        bool isInput = inputSet.find(i) != inputSet.end();
        bool isOutput = outputSet.find(i) != outputSet.end();
        bool isBias = false;
        
        if (!isInput && !isOutput) {
            bool hasOutgoingConnections = false;
            bool hasIncomingConnections = false;
            
            for (const auto& conn : connections) {
                if (conn._connectionGeneAttribute.enabled) {
                    if (conn._sourceNodeIndex == i) hasOutgoingConnections = true;
                    if (conn._targetNodeIndex == i) hasIncomingConnections = true;
                }
            }
            
            isBias = hasOutgoingConnections && !hasIncomingConnections;
        }
        
        // Skip bias nodes - they keep their value of 1.0
        if (isBias) {
            if (detailedDebug) {
                LOG_DEBUG("Skipping bias node[{}] (keeping value {:.3f})", i, nodeValues[i]);
            }
            continue;
        }
        
        // Apply activation function to hidden and output nodes
        double oldValue = nodeValues[i];
        if (nodes[i].activationType == ActivationType::SIGMOID) {
            nodeValues[i] = 1.0 / (1.0 + std::exp(-nodeInputs[i]));
        } else {
            nodeValues[i] = nodeInputs[i]; // No activation
        }
        
        if (detailedDebug) {
            LOG_DEBUG("Activate node[{}]: input={:.3f} -> output={:.3f} (was {:.3f})", 
                     i, nodeInputs[i], nodeValues[i], oldValue);
        }
    }
    
    // Return output node value
    if (!phenotype._outputIndices.empty()) {
        size_t outputIdx = phenotype._outputIndices[0];
        if (detailedDebug) {
            LOG_DEBUG("Final output: node[{}] = {:.3f}", outputIdx, nodeValues[outputIdx]);
            detailedDebug = false; // Only debug first network
        }
        return nodeValues[outputIdx];
    }
    
    return 0.0;
}

int main(int argc, char* argv[])
{
    configureLoggers(argc, argv);
    
    LOG_INFO("Starting NEAT XOR Evolution");
    
    try {
        // Evolution parameters
        const uint32_t populationSize = 500;
        const uint32_t maxGenerations = 40;
        const uint32_t randomSeed = 12345;
        
        // Create input node attributes (2 inputs for XOR)
        std::vector<NodeGeneAttributes> inputAttributes = {
            {ActivationType::NONE},  // Input 1
            {ActivationType::NONE}   // Input 2
        };
        
        // Create output node attributes (1 output for XOR)
        std::vector<NodeGeneAttributes> outputAttributes = {
            {ActivationType::SIGMOID}  // Output with sigmoid activation
        };
        
        // Create bias attributes (bias->output connection for XOR)
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttributes;
        biasAttributes[0] = {0.0, true}; // Bias to output with initial weight 0.0
        
        // Create initialization parameters
        Operator::InitParams initParams(
            inputAttributes,
            outputAttributes,
            biasAttributes,
            Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS
        );
        
        // Create plot elites parameters
        Population::PlotElitesParams eliteParams(
            0.2,  // elitePercentage - 20% of each species
            1,    // minimumElitesPerSpecies
            5     // maximumElitesPerSpecies
        );
        
        // Create plot crossover parameters
        Population::PlotCrossoverParams crossoverParams(
            0.5,  // topPerformerPercentage - top 50% can be parents
            2,    // baseCrossoversPerSpecies
            0.3,  // crossoverScalingFactor
            2     // minimumSpeciesSize for crossover
        );
        
        // Create crossover operator parameters
        Operator::CrossoverParams crossoverOperatorParams(
            0.25  // disabledGeneReactivationProbability
        );
        
        // Create cycle detection parameters
        Operator::CycleDetectionParams cycleDetectionParams;
        
        // Create node mutation parameters
        NodeGeneAttributes nodeAttribs{ActivationType::SIGMOID};
        Operator::NodeMutationParams nodeMutationParams(nodeAttribs);
        
        // Create connection mutation parameters
        Operator::ConnectionMutationParams connectionMutationParams(
            1.0,  // weight - higher initial weight for stronger signals
            4.0,  // weightRange - wider range for more exploration
            Operator::ConnectionMutationParams::NetworkTopology::FEED_FORWARD
        );
        
        // Create connection reactivation parameters
        Operator::ConnectionReactivationParams connectionReactivationParams(
            Operator::ConnectionReactivationParams::SelectionStrategy::RANDOM
        );
        
        // Create dynamic data update parameters
        Population::DynamicDataUpdateParams updateParams(
            5, // maxGenomePendingEliminationLimit
            20, // maxSpeciesPendingEliminationRating
            0.5, // speciesElitePlacementProtectionPercentage
            0.5, // speciesPendingEliminationPercentage
            0.3, // genomesPendingEliminationPercentage
            5, // equilibriumSpeciesCount
            populationSize // targetPopulationSize
        );
        
        // Create compatibility distance parameters
        Operator::CompatibilityDistanceParams compatibilityParams(
            1.0,  // c1 - excess gene coefficient
            1.0,  // c2 - disjoint gene coefficient
            0.4,  // c3 - weight difference coefficient
            8.0  // threshold - compatibility threshold (drastically increased to prevent species explosion)
        );
        
        // Create fitness strategy
        Analysis::SingleSimpleFitnessParams<DoubleFitnessResult> fitnessParams;
        fitnessParams.evaluationFunction = evaluateXORFitness;
        auto fitnessStrategy = std::make_unique<Analysis::SingleSimpleFitnessStrategy<DoubleFitnessResult>>(fitnessParams);
        
        // Create repair operator parameters (using defaults)
        Operator::RepairOperatorParams repairParams(2);
        
        // Create weight mutation parameters for better exploration
        Operator::WeightMutationParams weightMutationParams(
            0.8,  // perturbationRate - still high but allow more replacements
            0.2,  // replacementRate - increased for more dramatic changes
            1.0,  // perturbationStrength - stronger perturbations
            4.0,  // weightRange - wider range matching connection mutations
            Operator::WeightMutationParams::MutationType::MIXED
        );
        
        // Create mutation probability parameters optimized for XOR complexity
        Evolution::MutationProbabilityParams mutationParams(
            0.70,  // Weight mutation probability (70% - reduced to allow more structural mutations)
            0.10,  // Node mutation probability (10% - significantly increased for hidden nodes)
            0.15,  // Connection mutation probability (15% - increased for new connections)
            0.05   // Connection reactivation probability (5% - much higher to enable disabled connections)
        );
        
        // Create evolution prototype
        Evolution::EvolutionPrototype<DoubleFitnessResult> evolution(
            std::move(fitnessStrategy),
            populationSize,
            initParams,
            updateParams,
            eliteParams,
            crossoverParams,
            compatibilityParams,
            repairParams,
            mutationParams,
            weightMutationParams,
            crossoverOperatorParams,
            cycleDetectionParams,
            nodeMutationParams,
            connectionMutationParams,
            connectionReactivationParams,
            randomSeed
        );
        
        LOG_INFO("Starting evolution with {} individuals for {} generations", 
                 populationSize, maxGenerations);
        
        // Run evolution
        auto start = std::chrono::high_resolution_clock::now();
        auto results = evolution.run(maxGenerations);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        // Display results
        LOG_INFO("Evolution completed in {} ms", duration.count());
        LOG_INFO("Best fitness achieved: {}", results.getBestFitness().getValue());
        LOG_INFO("Final population size: {}", results.getFinalPopulation().size());
        
        // Display the best solution
        auto bestGenome = results.getBestGenome();
        
        std::cout << "\n=== XOR Problem Solution ===\n";
        std::cout << "Best fitness: " << results.getBestFitness().getValue() << "\n\n";
        
        // Display genome structure
        std::stringstream genomeOutput;
        Operator::displayGenome(bestGenome, genomeOutput);
        std::cout << genomeOutput.str();
        
        // Display phenotype structure
        std::stringstream phenotypeOutput;
        Operator::displayPhenotype(bestGenome, phenotypeOutput);
        std::cout << phenotypeOutput.str();
        
        LOG_INFO("NEAT XOR Evolution completed successfully");
        
    } catch (const std::exception& e) {
        LOG_ERROR("Evolution failed with exception: {}", e.what());
        return 1;
    }

    return 0;
}
