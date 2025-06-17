
#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>

#include "logger/Logger.hpp"

// Core NEAT components
#include "data/Genome.hpp"
#include "data/NodeGene.hpp"
#include "data/ConnectionGene.hpp"
#include "operator/Init.hpp"
#include "operator/DisplayGenome.hpp"
#include "operator/DisplayPhenotype.hpp"
#include "population/GenerationPlannerParams.hpp"
#include "population/DynamicDataUpdate.hpp"
#include "operator/CompatibilityDistance.hpp"
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
    
    // Test each XOR input pattern
    for (size_t i = 0; i < XOR_INPUTS.size(); ++i) {
        const auto& inputs = XOR_INPUTS[i];
        double expectedOutput = XOR_OUTPUTS[i];
        
        // TODO: Once network activation is available, replace this simulation
        // For now, simulate network behavior based on structure
        double networkOutput = simulateNetworkOutput(phenotype, inputs);
        
        // Calculate squared error
        double error = networkOutput - expectedOutput;
        totalError += error * error;
    }
    
    // Convert error to fitness (lower error = higher fitness)
    double performanceFitness = 4.0 - totalError; // Max fitness is 4.0 (perfect)
    
    // Add complexity penalization
    size_t nodeCount = phenotype._nodeGeneAttributes.size();
    size_t connectionCount = 0;
    
    // Count only enabled connections
    for (const auto& conn : phenotype._orderedConnections) {
        if (conn._connectionGeneAttribute.enabled) {
            connectionCount++;
        }
    }
    
    // Penalize excessive complexity (minimal network: 4 nodes, 3 connections)
    double complexityPenalty = (nodeCount > 4 ? (nodeCount - 4) * 0.1 : 0.0) + 
                              (connectionCount > 3 ? (connectionCount - 3) * 0.05 : 0.0);
    
    fitness = std::max(0.0, performanceFitness - complexityPenalty); // Ensure non-negative fitness
    
    return DoubleFitnessResult(fitness);
}

// Simulate network output based on structure (placeholder until network activation is implemented)
double simulateNetworkOutput(const Genome::Phenotype& phenotype, const std::vector<double>& inputs) {
    // Simple simulation: just return a weighted sum of inputs through direct connections
    // This is a placeholder - replace with proper network activation when available
    
    double output = 0.0;
    const auto& connections = phenotype._orderedConnections;
    
    // Find direct connections from inputs to output
    for (const auto& conn : connections) {
        if (conn._connectionGeneAttribute.enabled) {
            // Check if this is an input->output connection
            if (conn._sourceNodeIndex < inputs.size() && conn._targetNodeIndex == phenotype._outputIndices[0]) {
                output += inputs[conn._sourceNodeIndex] * conn._connectionGeneAttribute.weight;
            }
            // Add bias if source is bias node (assuming bias node index is 2)
            else if (conn._sourceNodeIndex == 2 && conn._targetNodeIndex == phenotype._outputIndices[0]) {
                output += 1.0 * conn._connectionGeneAttribute.weight;
            }
        }
    }
    
    // Apply sigmoid activation
    output = 1.0 / (1.0 + std::exp(-output));
    
    return output;
}

int main(int argc, char* argv[])
{
    configureLoggers(argc, argv);
    
    LOG_INFO("Starting NEAT XOR Evolution");
    
    try {
        // Evolution parameters
        const uint32_t populationSize = 200;
        const uint32_t maxGenerations = 20;
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
        
        // Create bias attributes (empty for simple XOR)
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttributes;
        
        // Create initialization parameters
        Operator::InitParams initParams(
            inputAttributes,
            outputAttributes,
            biasAttributes,
            Operator::InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS
        );
        
        // Create generation planner parameters
        auto plannerParams = Population::GenerationPlannerParamsFactory::createBalanced(
            populationSize, randomSeed
        );
        
        // Create dynamic data update parameters
        Population::DynamicDataUpdateParams updateParams(
            5,    // maxProtectionLimit
            3,    // maxSpeciesProtectionRating  
            0.3,  // protectedTierPercentage (30%)
            1     // worstSpeciesCount
        );
        
        // Create compatibility distance parameters
        Operator::CompatibilityDistanceParams compatibilityParams(
            1.0,  // c1 - excess gene coefficient
            1.0,  // c2 - disjoint gene coefficient
            0.4,  // c3 - weight difference coefficient
            3.0   // threshold - compatibility threshold
        );
        
        // Create fitness strategy
        Analysis::SingleSimpleFitnessParams<DoubleFitnessResult> fitnessParams;
        fitnessParams.evaluationFunction = evaluateXORFitness;
        auto fitnessStrategy = std::make_unique<Analysis::SingleSimpleFitnessStrategy<DoubleFitnessResult>>(fitnessParams);
        
        // Create evolution prototype
        Evolution::EvolutionPrototype<DoubleFitnessResult> evolution(
            std::move(fitnessStrategy),
            populationSize,
            initParams,
            plannerParams,
            updateParams,
            compatibilityParams,
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
