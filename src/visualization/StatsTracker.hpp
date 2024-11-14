// stats_tracker.hpp
#pragma once

#include <string>

#include "core/Genome.hpp"

namespace neat {
namespace core {
class Population;
}
namespace visualization {

class StatsTracker {
public:
    struct GenerationStats {
        int32_t generation;
        double bestFitness;
        double averageFitness;
        double worstFitness;
        int32_t speciesCount;
        int32_t totalNodes;
        int32_t totalConnections;
        std::vector<double> speciesFitnesses;
        std::vector<size_t> speciesSizes;
    };

    void recordGeneration(const core::Population& population);
    void recordBestGenome(const core::Genome& genome);
    
    // Generate various plots
    std::string generateFitnessPlot();
    std::string generateSpeciesPlot();
    std::string generateComplexityPlot();
    
    // Export data
    void exportToCSV(const std::string& filename);
    void saveAllPlots(const std::string& directory);
    
    const std::vector<GenerationStats>& getStats() const { return generationStats; }
    const std::vector<core::Genome>& getBestGenomes() const { return bestGenomes; }

private:
    std::string generatePlot(const std::string& title,
                            const std::string& xLabel,
                            const std::string& yLabel,
                            const std::vector<double>& data);

    std::vector<GenerationStats> generationStats;
    std::vector<core::Genome> bestGenomes;
};

}
}
