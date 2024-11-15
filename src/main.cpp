#include "core/NEAT.hpp"
#include "visualization/NetworkRenderer.hpp"
#include "core/Gene.hpp"
#include "core/Species.hpp"
#include <iostream>
#include <vector>

using namespace neat;

// XOR truth table inputs/outputs
const std::vector<std::vector<double>> XOR_INPUTS = {
    {0, 0}, {0, 1}, {1, 0}, {1, 1}
};
const std::vector<double> XOR_OUTPUTS = {0, 1, 1, 0};

double evaluateXORFitness(const core::Genome& genome) {
    double fitness = 4.0;  // Start with max fitness and subtract errors
    
    // Get number of inputs the network expects
    size_t inputCount = std::count_if(genome.getNodes().begin(), genome.getNodes().end(),
        [](const auto& pair) { return pair.second == neat::core::ENodeType::INPUT; });
    
    std::vector<double> inputs(inputCount);  // Will include bias if present
    
    for (size_t i = 0; i < XOR_INPUTS.size(); i++) {
        // Copy the XOR inputs
        for (size_t j = 0; j < XOR_INPUTS[i].size(); j++) {
            inputs[j] = XOR_INPUTS[i][j];
        }
        
        // If we have an extra input (bias), set it to 1.0
        if (inputs.size() > XOR_INPUTS[i].size()) {
            inputs.back() = 1.0;
        }
        
        auto outputs = genome.activate(inputs);
        if (outputs.empty()) {
            return 0.0;  // Handle error case
        }
        
        double error = std::abs(outputs[0] - XOR_OUTPUTS[i]);
        fitness -= error * error;  // Square error
    }
    
    return std::max(0.0, fitness);
}

void generationCallback(int32_t generation, const core::Population& population) {
    const auto& best = population.getBestGenome();
    std::cout << "Generation " << generation 
              << " Best Fitness: " << best.getFitness() << std::endl;
}

int main() {
    core::NEAT::Config config(2, 1);
    config.populationConfig.populationSize = 150;
    
    // Print configuration
    std::cout << "Creating NEAT with configuration:" << std::endl
              << "Inputs: " << config.inputSize << std::endl
              << "Outputs: " << config.outputSize << std::endl
              << "Population size: " << config.populationConfig.populationSize << std::endl;
    
    core::NEAT neat(config);
    
    // Run evolution
    neat.evolve(100, evaluateXORFitness, generationCallback);
    
    // Visualize best solution
    visualization::NetworkRenderer::Config rendererConfig;
	visualization::NetworkRenderer renderer(rendererConfig);
    auto svg = renderer.renderToSVG(neat.getBestGenome());
    std::cout << "\nBest network visualization:\n" << svg << std::endl;
    
    return 0;
}
