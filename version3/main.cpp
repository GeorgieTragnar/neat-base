
#include <fstream>
#include <iostream>
#include <filesystem>
#include <map>

#include "logger/Logger.hpp"

#include "data/NodeGene.hpp"
static auto logger = LOGGER("main");

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

// double evaluateXORFitness(const core::Genome& genome) {
//     double fitness = 4.0;  // Start with max fitness and subtract errors
    
//     // Get number of inputs the network expects
//     size_t inputCount = std::count_if(genome.getNodes().begin(), genome.getNodes().end(),
//         [](const auto& pair) { return pair.second == neat::core::ENodeType::INPUT; });
    
//     std::vector<double> inputs(inputCount);  // Will include bias if present
    
//     for (size_t i = 0; i < XOR_INPUTS.size(); i++) {
//         // Copy the XOR inputs
//         for (size_t j = 0; j < XOR_INPUTS[i].size(); j++) {
//             inputs[j] = XOR_INPUTS[i][j];
//         }
        
//         // If we have an extra input (bias), set it to 1.0
//         if (inputs.size() > XOR_INPUTS[i].size()) {
//             inputs.back() = 1.0;
//         }
        
//         auto outputs = genome.activate(inputs);
//         if (outputs.empty()) {
//             return 0.0;  // Handle error case
//         }
        
//         double error = std::abs(outputs[0] - XOR_OUTPUTS[i]);
//         fitness -= error * error;  // Square error
//     }
    
//     return std::max(0.0, fitness);
// }

int main(int argc, char* argv[])
{
    configureLoggers(argc, argv);


}
