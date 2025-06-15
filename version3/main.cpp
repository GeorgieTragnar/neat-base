
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
#include "population/GenerationPlannerParams.hpp"
#include "population/DynamicDataUpdate.hpp"
#include "operator/CompatibilityDistance.hpp"
#include "analysis/strategies/SingleSimpleFitnessStrategy.hpp"

// Evolution prototype
#include "evolution/EvolutionPrototype.hpp"

static auto logger = LOGGER("main");

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

// XOR fitness evaluation function
DoubleFitnessResult evaluateXORFitness(const Genome::Phenotype& phenotype) {
    // For now, return a simple fitness based on network complexity
    // TODO: Implement proper network activation when available
    double fitness = 1.0;
    
    // Simple fitness based on network structure
    size_t nodeCount = phenotype._nodeGeneAttributes.size();
    size_t connectionCount = phenotype._orderedConnections.size();
    
    // Reward having some connections and nodes
    fitness += connectionCount * 0.1;
    fitness += nodeCount * 0.05;
    
    return DoubleFitnessResult(fitness);
}

int main(int argc, char* argv[])
{
    configureLoggers(argc, argv);
    
    LOG_INFO("Starting NEAT XOR Evolution");
    
    try {
        // Evolution parameters
        const uint32_t populationSize = 150;
        const uint32_t maxGenerations = 300;
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
        
        // Test the best solution
        auto bestGenome = results.getBestGenome();
        // const auto& bestPhenotype = bestGenome.get_phenotype();
        
        // std::cout << "\n=== XOR Problem Solution ===\n";
        // std::cout << "Best fitness: " << results.getBestFitness() << "\n";
        // std::cout << "Testing best network:\n";
        
        // // Show network structure instead of activation results
        // std::cout << "Network structure:\n";
        // std::cout << "Nodes: " << bestPhenotype._nodeGeneAttributes.size() << "\n";
        // std::cout << "Connections: " << bestPhenotype._orderedConnections.size() << "\n";
        // std::cout << "Inputs: " << bestPhenotype._inputIndices.size() << "\n";
        // std::cout << "Outputs: " << bestPhenotype._outputIndices.size() << "\n";
        
        LOG_INFO("NEAT XOR Evolution completed successfully");
        
    } catch (const std::exception& e) {
        LOG_ERROR("Evolution failed with exception: {}", e.what());
        return 1;
    }

    return 0;
}
