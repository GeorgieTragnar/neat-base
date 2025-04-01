#include "core/NEAT.hpp"
#include "visualization/NetworkRenderer.hpp"
#include "core/Gene.hpp"
#include "core/Species.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>

#include "logger/Logger.hpp"
static auto logger = LOGGER("main");

using namespace neat;


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

// Proper sine function inputs/outputs
const std::vector<std::vector<double>> SIN_INPUTS = {
    {0.0}, {0.2}, {0.4}, {0.6}, {0.8}, {1.0}, {1.2}, {1.4}, {1.6}, {1.8}, {2.0},
    {2.2}, {2.4}, {2.6}, {2.8}, {3.0}, {3.2}, {3.4}, {3.6}, {3.8}, {4.0}, {4.2},
    {4.4}, {4.6}, {4.8}, {5.0}, {5.2}, {5.4}, {5.6}, {5.8}, {6.0}, {6.2}
};

// Precalculated sin values - no need to compute at runtime
const std::vector<double> SIN_OUTPUTS = {
    0.0, 0.199, 0.389, 0.565, 0.717, 0.841, 0.932, 0.985, 0.999, 0.973, 0.909,
    0.808, 0.675, 0.515, 0.335, 0.141, -0.058, -0.256, -0.443, -0.616, -0.757,
    -0.871, -0.951, -0.994, -0.997, -0.959, -0.883, -0.772, -0.631, -0.464, -0.279, -0.083
};

double evaluateSinFitness(core::Genome& genome) {
    // Validate network dimensions
    if (genome.getConfig().inputSize != 1 || genome.getConfig().outputSize != 1) {
        throw std::runtime_error("Network must have 1 input and 1 output for sin function");
    }
    
    double maxError = SIN_INPUTS.size() * 2.0; // Max error of 2.0 per point (-1 to 1 range)
    double totalError = 0.0;
    
    std::vector<double> inputs(1);
    
    for (size_t i = 0; i < SIN_INPUTS.size(); i++) {
        // Scale inputs to a reasonable range for the network
        // Most activation functions work well in [-1,1] range
        inputs[0] = (SIN_INPUTS[i][0] / 3.0) - 1.0; // Scale to [-1,1] range
        
        auto outputs = genome.activate(inputs);
        if (outputs.empty()) return 0.0;
        
        // Calculate absolute error
        double error = std::abs(outputs[0] - SIN_OUTPUTS[i]);
        totalError += error;
    }
    
    // Convert error to fitness (0 to 1 range)
    double fitness = 1.0 - (totalError / maxError);
    
    return std::max(0.0, fitness);
}

// XOR truth table inputs/outputs
const std::vector<std::vector<double>> XOR_INPUTS = {
    {0, 0}, {0, 1}, {1, 0}, {1, 1}
};
const std::vector<double> XOR_OUTPUTS = {0, 1, 1, 0};

// double evaluateXORFitness(core::Genome& genome) {
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

void generationCallback(int32_t generation, const core::Population& population) {
    const auto& best = population.getBestGenome();
    LOG_INFO("Generation {} Best Fitness {}", generation, best.getFitness());
}

int main(int argc, char* argv[]) {
    // Configure loggers from command line arguments
    configureLoggers(argc, argv);

    core::NEAT::Config config(1, 1);
    
    LOG_INFO("Creating NEAT with configuration:");
    LOG_INFO("Inputs: {}", config.inputSize);
    LOG_INFO("Outputs: {}", config.outputSize);
    
    core::NEAT neat(config);
    
    // Create outputs directory if it doesn't exist
    std::filesystem::create_directories("outputs");
    
    // Get current time for folder name
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "outputs/run_" << std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S");
    std::string outputDir = ss.str();
    std::filesystem::create_directories(outputDir);
    
    // Run evolution
    neat.evolve(50, evaluateSinFitness, generationCallback);
    
    // Get the best genome
    const auto& bestGenome = neat.getBestGenome();

    // Visualize best solution
    visualization::NetworkRenderer::Config rendererConfig;
	visualization::NetworkRenderer renderer(rendererConfig);
    auto svg = renderer.renderToSVG(bestGenome);
	
            
            // Save SVG file
            std::ofstream svgFile(outputDir + "/network.svg");
            svgFile << svg;
            svgFile.close();
            
            // Create HTML viewer
            std::ofstream htmlFile(outputDir + "/index.html");
            // In main.cpp where you generate the HTML
            htmlFile << R"(<!DOCTYPE html>
<html>
<head>
    <title>NEAT Network Visualization</title>
    <style>
        body {
            display: flex;
            flex-direction: column;
            align-items: center;
            min-height: 100vh;
            margin: 0;
            background: #f0f0f0;
            font-family: Arial, sans-serif;
        }
        .header {
            width: 100%;
            background: #2196F3;
            color: white;
            padding: 1rem;
            text-align: center;
            margin-bottom: 2rem;
        }
        .svg-container {
            background: white;
            padding: 20px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            margin-bottom: 2rem;
        }
        .info {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
            margin-bottom: 2rem;
            max-width: 800px;
            width: 100%;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1rem;
        }
        .stats-item {
            background: #f8f9fa;
            padding: 1rem;
            border-radius: 4px;
        }
        .stats-item h3 {
            margin-top: 0;
            color: #2196F3;
            border-bottom: 2px solid #2196F3;
            padding-bottom: 0.5rem;
        }
        .connection-list {
            margin-top: 1rem;
            font-family: monospace;
            background: #f8f9fa;
            padding: 1rem;
            border-radius: 4px;
        }
        .stat-value {
            font-weight: bold;
            color: #2196F3;
        }
        .connection-detail {
            display: grid;
            grid-template-columns: auto 1fr;
            gap: 0.5rem;
            align-items: center;
            padding: 0.25rem 0;
        }
        .connection-weight {
            color: #2196F3;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>NEAT Network Visualization</h1>
        <p>Run completed at )" << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S") << R"(</p>
    </div>
    <div class="info">
        <h2>Run Information</h2>
        <div class="stats-grid">
            <div class="stats-item">
                <h3>Performance</h3>
                <p>Best Fitness: <span class="stat-value">)" << bestGenome.getFitness() << R"(</span></p>
                <p>Total Generations: <span class="stat-value">100</span></p>
            </div>
            <div class="stats-item">
                <h3>Network Structure</h3>
                <p>Input Nodes: <span class="stat-value">)" << 
                std::count_if(bestGenome.getNodes().begin(), bestGenome.getNodes().end(),
                    [](const auto& pair) { return pair.second == core::ENodeType::INPUT; }) << R"(</span></p>
                <p>Output Nodes: <span class="stat-value">)" << 
                std::count_if(bestGenome.getNodes().begin(), bestGenome.getNodes().end(),
                    [](const auto& pair) { return pair.second == core::ENodeType::OUTPUT; }) << R"(</span></p>
                <p>Hidden Nodes: <span class="stat-value">)" << 
                std::count_if(bestGenome.getNodes().begin(), bestGenome.getNodes().end(),
                    [](const auto& pair) { return pair.second == core::ENodeType::HIDDEN; }) << R"(</span></p>
                <p>Enabled Connections: <span class="stat-value">)" << 
                std::count_if(bestGenome.getGenes().begin(), bestGenome.getGenes().end(),
                    [](const auto& gene) { return gene.enabled; }) << R"(</span></p>
                <p>Total Connections: <span class="stat-value">)" << bestGenome.getGenes().size() << R"(</span></p>
            </div>
        </div>
        
        <div class="connection-list">
            <h3>Connection Details</h3>
            <div style="max-height: 200px; overflow-y: auto;">
            <pre>)";
            
            // Add connection details
            for (const auto& gene : bestGenome.getGenes()) {
                htmlFile << (gene.enabled ? "✓" : "✗") << " ";
                htmlFile << "Node " << std::setw(3) << gene.inputNode << " → Node " << std::setw(3) << gene.outputNode 
                        << " (weight: " << std::fixed << std::setprecision(2) 
                        << std::setw(8) << gene.weight << ")" 
                        << (gene.enabled ? "" : " [disabled]") << "\n";
            }
            
            htmlFile << R"(</pre>
            </div>
        </div>

        <div class="connection-list">
            <h3>Network Analysis</h3>
            <div style="max-height: 200px; overflow-y: auto;">
            <pre>)";

            // Calculate connected nodes
            std::set<int32_t> connectedNodes;
            for (const auto& gene : bestGenome.getGenes()) {
                if (gene.enabled) {
                    connectedNodes.insert(gene.inputNode);
                    connectedNodes.insert(gene.outputNode);
                }
            }

            // Print connected nodes
            htmlFile << "Connected nodes:\n";
            for (int32_t id : connectedNodes) {
                const auto& type = bestGenome.getNodes().at(id);
                htmlFile << "  Node " << std::setw(3) << id << ": "
                        << (type == core::ENodeType::INPUT ? "INPUT" :
                            type == core::ENodeType::OUTPUT ? "OUTPUT" :
                            type == core::ENodeType::BIAS ? "BIAS" : "HIDDEN")
                        << "\n";
            }

            // Print disconnected nodes
            htmlFile << "\nDisconnected nodes:\n";
            for (const auto& [id, type] : bestGenome.getNodes()) {
                if (!connectedNodes.count(id)) {
                    htmlFile << "  Node " << std::setw(3) << id << ": "
                            << (type == core::ENodeType::INPUT ? "INPUT" :
                                type == core::ENodeType::OUTPUT ? "OUTPUT" :
                                type == core::ENodeType::BIAS ? "BIAS" : "HIDDEN")
                            << "\n";
                }
            }

            htmlFile << R"(</pre>
            </div>
        </div>
    </div>
    <div class="svg-container">)"
        << svg <<
    R"(</div>
</body>
</html>)";
            htmlFile.close();
            
            // Create bash script to open visualization
            std::ofstream scriptFile(outputDir + "/view.sh");
            scriptFile << "#!/bin/bash\n"
                      << "xdg-open \"$(dirname \"$0\")/index.html\"\n";
            scriptFile.close();
            
            // Make script executable
            std::filesystem::permissions(
                outputDir + "/view.sh",
                std::filesystem::perms::owner_exec | 
                std::filesystem::perms::group_exec |
                std::filesystem::perms::others_exec,
                std::filesystem::perm_options::add
            );
            
            LOG_INFO("Visualization saved to {}", outputDir);
            LOG_INFO("Run './view.sh' in that directory to view the results");
        
    
    return 0;
}
