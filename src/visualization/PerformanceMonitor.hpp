// performance_monitor.hpp
#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>

namespace neat {
namespace visualization {

class PerformanceMonitor {
public:
    struct Metrics {
        double evaluationTime;
        double mutationTime;
        double crossoverTime;
        double speciationTime;
        double totalTime;
        size_t activeSpecies;
        size_t totalEvaluations;
        size_t innovations;
        double memoryUsage;
    };

    void startGeneration();
    void endGeneration();
    
    void startPhase(const std::string& phase);
    void endPhase(const std::string& phase);
    
    void recordMetric(const std::string& name, double value);
    
    std::string generateReport();
    void exportMetrics(const std::string& filename);
    
    const std::vector<Metrics>& getMetrics() const { return generationMetrics; }

private:
    std::vector<Metrics> generationMetrics;
    std::unordered_map<std::string, std::chrono::steady_clock::time_point> phaseTimers;
    Metrics currentMetrics;
};

}
}
