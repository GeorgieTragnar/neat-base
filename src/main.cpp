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
    std::vector<double> inputs(3);  // 2 inputs + bias
    inputs[2] = 1.0;  // Bias input
    
    for (size_t i = 0; i < XOR_INPUTS.size(); i++) {
        inputs[0] = XOR_INPUTS[i][0];
        inputs[1] = XOR_INPUTS[i][1];
        
        auto outputs = genome.activate(inputs);
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
    core::NEAT::Config config;
    config.inputSize = 3;  // 2 inputs + bias
    config.outputSize = 1;
    config.populationConfig.populationSize = 150;
    
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
