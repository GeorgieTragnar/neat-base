#include "StatsTracker.hpp"
#include "core/Population.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>

namespace neat {
namespace visualization {

void StatsTracker::recordGeneration(const core::Population& population) {
    GenerationStats stats;
    stats.generation = generationStats.size();
    
    // Calculate fitness statistics
    stats.bestFitness = std::numeric_limits<double>::lowest();
    stats.worstFitness = std::numeric_limits<double>::max();
    double sumFitness = 0.0;
    
    for (size_t i = 0; i < population.size(); ++i) {
        const auto& genome = population.getGenome(i);
        double fitness = genome.getFitness();
        stats.bestFitness = std::max(stats.bestFitness, fitness);
        stats.worstFitness = std::min(stats.worstFitness, fitness);
        sumFitness += fitness;
    }
    
    stats.averageFitness = sumFitness / population.size();
    
    // Species statistics
    const auto& species = population.getSpecies();
    stats.speciesCount = species.size();
    
    for (const auto& s : species) {
        stats.speciesSizes.push_back(s.getMembers().size());
        double speciesBestFitness = 0.0;
        for (auto idx : s.getMembers()) {
            speciesBestFitness = std::max(speciesBestFitness, 
                population.getGenome(idx).getFitness());
        }
        stats.speciesFitnesses.push_back(speciesBestFitness);
    }
    
    // Complexity statistics
    stats.totalNodes = 0;
    stats.totalConnections = 0;
    for (size_t i = 0; i < population.size(); ++i) {
        const auto& genome = population.getGenome(i);
        stats.totalNodes += genome.getNodes().size();
        stats.totalConnections += genome.getGenes().size();
    }
    
    generationStats.push_back(stats);
}

void StatsTracker::recordBestGenome(const core::Genome& genome) {
    bestGenomes.push_back(genome);
}

std::string StatsTracker::generatePlot(const std::string& title,
                                     const std::string& xLabel,
                                     const std::string& yLabel,
                                     const std::vector<double>& data) {
    if (data.empty()) return "";
    
    std::ostringstream ss;
    ss << "<svg width=\"800\" height=\"400\">\n";
    
    // Calculate scaling
    double maxY = *std::max_element(data.begin(), data.end());
    double minY = *std::min_element(data.begin(), data.end());
    double range = maxY - minY;
    double scale = 350.0 / range;
    
    // Draw axes
    ss << "<line x1=\"50\" y1=\"350\" x2=\"750\" y2=\"350\" stroke=\"black\" />\n";
    ss << "<line x1=\"50\" y1=\"350\" x2=\"50\" y2=\"50\" stroke=\"black\" />\n";
    
    // Draw labels
    ss << "<text x=\"400\" y=\"390\" text-anchor=\"middle\">" << xLabel << "</text>\n";
    ss << "<text x=\"20\" y=\"200\" transform=\"rotate(-90, 20, 200)\">" << yLabel << "</text>\n";
    ss << "<text x=\"400\" y=\"30\" text-anchor=\"middle\">" << title << "</text>\n";
    
    // Draw data points and lines
    ss << "<path d=\"M";
    for (size_t i = 0; i < data.size(); ++i) {
        double x = 50 + (700.0 * i / (data.size() - 1));
        double y = 350 - (data[i] - minY) * scale;
        ss << x << "," << y << " ";
    }
    ss << "\" fill=\"none\" stroke=\"blue\" stroke-width=\"2\" />\n";
    
    ss << "</svg>";
    return ss.str();
}

std::string StatsTracker::generateFitnessPlot() {
    std::vector<double> bestFitnesses;
    for (const auto& stats : generationStats) {
        bestFitnesses.push_back(stats.bestFitness);
    }
    return generatePlot("Fitness Over Time", "Generation", "Fitness", bestFitnesses);
}

std::string StatsTracker::generateSpeciesPlot() {
    std::vector<double> speciesCounts;
    for (const auto& stats : generationStats) {
        speciesCounts.push_back(stats.speciesCount);
    }
    return generatePlot("Species Count Over Time", "Generation", "Species", speciesCounts);
}

std::string StatsTracker::generateComplexityPlot() {
    std::vector<double> complexity;
    for (const auto& stats : generationStats) {
        complexity.push_back(stats.totalNodes + stats.totalConnections);
    }
    return generatePlot("Network Complexity", "Generation", "Total Elements", complexity);
}

void StatsTracker::exportToCSV(const std::string& filename) {
    std::ofstream file(filename);
    if (!file) return;
    
    file << "Generation,BestFitness,AvgFitness,WorstFitness,Species,Nodes,Connections\n";
    
    for (const auto& stats : generationStats) {
        file << stats.generation << ","
             << stats.bestFitness << ","
             << stats.averageFitness << ","
             << stats.worstFitness << ","
             << stats.speciesCount << ","
             << stats.totalNodes << ","
             << stats.totalConnections << "\n";
    }
}

void StatsTracker::saveAllPlots(const std::string& directory) {
    std::filesystem::create_directories(directory);
    
    std::ofstream fitness(directory + "/fitness.svg");
    fitness << generateFitnessPlot();
    
    std::ofstream species(directory + "/species.svg");
    species << generateSpeciesPlot();
    
    std::ofstream complexity(directory + "/complexity.svg");
    complexity << generateComplexityPlot();
}

}
}
