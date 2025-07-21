
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
#include "population/Init.hpp"
#include "evolution/WeightMutation.hpp"
#include "display/DisplayGenome.hpp"
#include "display/DisplayPhenotype.hpp"
#include "population/DynamicDataUpdate.hpp"
#include "population/PlotElites.hpp"
#include "evolution/PlotCrossover.hpp"
#include "analysis/CompatibilityDistance.hpp"
#include "evolution/Crossover.hpp"
#include "validation/CycleDetection.hpp"
#include "evolution/NodeMutation.hpp"
#include "evolution/ConnectionMutation.hpp"
#include "evolution/ConnectionReactivation.hpp"
#include "evolution/RepairOperator.hpp"
#include "analysis/strategies/SingleSimpleFitnessStrategy.hpp"
#include "phenotype/NetworkExecution.hpp"

// Evolution prototype
#include "EvolutionPrototype.hpp"

static auto logger = LOGGER("evolution.EvolutionPrototype");

// User-defined fitness result type for XOR problem
class DoubleFitnessResult {
public:
    explicit DoubleFitnessResult(double penalizedValue = 0.0, double rawValue = 0.0) 
        : _penalizedValue(penalizedValue), _rawValue(rawValue) {}
    
    bool isBetterThan(const DoubleFitnessResult& other) const {
        return _penalizedValue > other._penalizedValue; // Higher is better, use penalized for comparison
    }
    
    bool isEqualTo(const DoubleFitnessResult& other) const {
        return std::abs(_penalizedValue - other._penalizedValue) < std::numeric_limits<double>::epsilon();
    }
    
    double getValue() const { return _penalizedValue; } // Return penalized for sorting/selection
    double getRawValue() const { return _rawValue; }   // Raw fitness without complexity penalty
    
    // Comparison operators for multimap ordering (use penalized value)
    bool operator<(const DoubleFitnessResult& other) const {
        return _penalizedValue < other._penalizedValue;
    }

private:
    double _penalizedValue; // Fitness with complexity penalty (used for selection)
    double _rawValue;       // Raw fitness without complexity penalty (for display)
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
double simulateNetworkOutput(const Phenotype& phenotype, const std::vector<double>& inputs);

// XOR fitness evaluation function
DoubleFitnessResult evaluateXORFitness(const Phenotype& phenotype) {
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
    
    // Debug output for high-fitness networks (higher threshold for success)
    if (performanceFitness > 3.9) {
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
    
    
    // Calculate raw fitness (without complexity penalty)
    double rawFitness;
    if (performanceFitness <= 3.0) {
        rawFitness = performanceFitness * 0.5; // Halve the fitness for poor performers
        fitness = rawFitness; // No complexity penalty for poor performers
    } else {
        // For structures showing signs of capability, quadruple the fitness
        double quadrupledFitness = performanceFitness * 4.0; // Max becomes 16
        rawFitness = quadrupledFitness; // Raw fitness without complexity penalty
        
        // Calculate complexity penalty with max of 12
        size_t nodeCount = phenotype._nodeGeneAttributes.size();
        size_t connectionCount = 0;
        
        // Count only enabled connections
        for (const auto& conn : phenotype._orderedConnections) {
            if (conn._connectionGeneAttribute.enabled) {
                connectionCount++;
            }
        }
        
        // Complexity penalty calculation (ideal: 4 nodes, 6 connections)
        double complexityPenalty = (nodeCount > 4 ? (nodeCount - 4) * 0.15 : 0.0) + 
                                  (connectionCount > 6 ? (connectionCount - 6) * 0.05 : 0.0);
        
        // Clamp complexity penalty to maximum of 12
        complexityPenalty = std::min(complexityPenalty, 12.0);
        
        // Apply complexity penalty for penalized fitness
        fitness = std::max(0.0, quadrupledFitness - complexityPenalty);
    }
    
    return DoubleFitnessResult(fitness, rawFitness);
}

// Use NetworkExecution operator for consistent network evaluation
double simulateNetworkOutput(const Phenotype& phenotype, const std::vector<double>& inputs) {
    // Create NetworkExecution parameters with debug output for high-fitness networks
    static bool debugOutput = true; // Enable debug for first few networks
    static int executionCount = 0;
    
    // Only debug first few executions or high-fitness networks
    bool shouldDebug = debugOutput && executionCount < 3;
    
    Operator::NetworkExecutionParams params(shouldDebug);
    
    // Execute network using the NetworkExecution operator
    auto outputs = Operator::networkExecution(phenotype, inputs, params);
    
    executionCount++;
    if (executionCount >= 3) {
        debugOutput = false; // Disable debug after first few networks
    }
    
    // Return first output (XOR has only one output)
    return outputs.empty() ? 0.0 : outputs[0];
}

int main(int argc, char* argv[])
{
    configureLoggers(argc, argv);
    
    LOG_INFO("Starting NEAT XOR Evolution");
    
    try {
        // Evolution parameters
        const uint32_t populationSize = 100;
        const uint32_t maxGenerations = 160;
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
        
        // Create bias attributes with random initial weight for better learning
        std::mt19937 rng(randomSeed);
        std::uniform_real_distribution<double> biasWeightDist(-1.0, 1.0);
        
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttributes;
        // biasAttributes[0] = {float(biasWeightDist(rng)), true}; // Random bias weight [-1.0, 1.0]
        
        // Create initialization parameters
        Operator::InitParams initParams(
            inputAttributes,
            outputAttributes,
            biasAttributes,
            Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS
        );
        
        // Create plot elites parameters
        Operator::PlotElitesParams eliteParams(
            0.01,  // elitePercentage - 20% of each species
            1,    // minimumElitesPerSpecies
            1     // maximumElitesPerSpecies
        );
        
        // Create plot crossover parameters
        Operator::PlotCrossoverParams crossoverParams(
            0.5,  // topPerformerPercentage - top 50% can be parents
            2,    // baseCrossoversPerSpecies
            0.05, // crossoverScalingFactor - reduced from 0.3 to 0.05 (5%)
            2     // minimumSpeciesSize for crossover
        );
        
        // Create crossover operator parameters
        Operator::CrossoverParams crossoverOperatorParams(
            0.25  // disabledGeneReactivationProbability
        );
        
        // Create cycle detection parameters
        Operator::CycleDetectionParams cycleDetectionParams(true, true);
        
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
        Operator::DynamicDataUpdateParams updateParams(
            2, // maxGenomePendingEliminationLimit
            10, // maxSpeciesPendingEliminationRating - reduced from 10 to 3 (faster elimination)
            1, // speciesElitePlacementProtectionPercentage
            0.9, // speciesPendingEliminationPercentage - increased from 0.8 to 0.95 (95%)
            0.6, // genomesPendingEliminationPercentage - increased from 0.3 to 0.7 (70%)
            5, // equilibriumSpeciesCount
            populationSize // targetPopulationSize
        );
        
        // Create compatibility distance parameters
        Operator::CompatibilityDistanceParams compatibilityParams(
            1.0,  // c1 - excess gene coefficient
            1.0,  // c2 - disjoint gene coefficient
            0.4,  // c3 - weight difference coefficient
            10.0  // threshold - compatibility threshold (drastically increased to prevent species explosion)
        );
        
        // Create fitness strategy
        Analysis::SingleSimpleFitnessParams<DoubleFitnessResult> fitnessParams;
        fitnessParams.evaluationFunction = evaluateXORFitness;
        auto fitnessStrategy = std::make_shared<Analysis::SingleSimpleFitnessStrategy<DoubleFitnessResult>>(fitnessParams);
        
        // Create repair operator parameters (using defaults)
        Operator::RepairOperatorParams repairParams(2);
        
        // Create weight mutation parameters for better exploration
        Operator::WeightMutationParams weightMutationParams(
            0.8,  // perturbationRate - still high but allow more replacements
            0.25, // replacementRate - increased for more dramatic changes
            1.5,  // perturbationStrength - stronger perturbations for better exploration
            4.0,  // weightRange - wider range matching connection mutations
            Operator::WeightMutationParams::MutationType::MIXED
        );
        
        // Create mutation probability parameters tuned for simpler solutions
        Evolution::MutationProbabilityParams mutationParams(
            0.80,  // Weight mutation probability (80% - increased for weight optimization)
            0.02,  // Node mutation probability (2% - further reduced to minimize structural bloat)
            0.18,  // Connection mutation probability (18% - increased for better connectivity exploration)
            0.05   // Connection reactivation probability (5% - same)
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
        std::cout << "Best fitness (with complexity penalty): " << results.getBestFitness().getValue() << "\n";
        std::cout << "Raw fitness (without complexity penalty): " << results.getBestFitness().getRawValue() << "\n\n";
        
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
