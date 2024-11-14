#include "PerformanceMonitor.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace neat {
namespace visualization {

void PerformanceMonitor::startGeneration() {
    currentMetrics = Metrics{};
    phaseTimers.clear();
    startPhase("total");
}

void PerformanceMonitor::endGeneration() {
    endPhase("total");
    currentMetrics.totalTime = currentMetrics.evaluationTime + 
                              currentMetrics.mutationTime + 
                              currentMetrics.crossoverTime + 
                              currentMetrics.speciationTime;
    generationMetrics.push_back(currentMetrics);
}

void PerformanceMonitor::startPhase(const std::string& phase) {
    phaseTimers[phase] = std::chrono::steady_clock::now();
}

void PerformanceMonitor::endPhase(const std::string& phase) {
    auto end = std::chrono::steady_clock::now();
    auto start = phaseTimers[phase];
    double duration = std::chrono::duration<double>(end - start).count();

    if (phase == "evaluation") currentMetrics.evaluationTime = duration;
    else if (phase == "mutation") currentMetrics.mutationTime = duration;
    else if (phase == "crossover") currentMetrics.crossoverTime = duration;
    else if (phase == "speciation") currentMetrics.speciationTime = duration;
}

void PerformanceMonitor::recordMetric(const std::string& name, double value) {
    if (name == "active_species") currentMetrics.activeSpecies = static_cast<size_t>(value);
    else if (name == "total_evaluations") currentMetrics.totalEvaluations = static_cast<size_t>(value);
    else if (name == "innovations") currentMetrics.innovations = static_cast<size_t>(value);
    else if (name == "memory_usage") currentMetrics.memoryUsage = value;
}

std::string PerformanceMonitor::generateReport() {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3);
    ss << "Performance Report\n";
    ss << "=================\n\n";

    if (generationMetrics.empty()) {
        ss << "No metrics recorded yet.\n";
        return ss.str();
    }

    // Calculate averages
    double avgEvalTime = 0, avgMutTime = 0, avgCrossTime = 0, avgSpecTime = 0;
    double avgMemUsage = 0;
    size_t avgSpecies = 0;

    for (const auto& metrics : generationMetrics) {
        avgEvalTime += metrics.evaluationTime;
        avgMutTime += metrics.mutationTime;
        avgCrossTime += metrics.crossoverTime;
        avgSpecTime += metrics.speciationTime;
        avgMemUsage += metrics.memoryUsage;
        avgSpecies += metrics.activeSpecies;
    }

    size_t count = generationMetrics.size();
    avgEvalTime /= count;
    avgMutTime /= count;
    avgCrossTime /= count;
    avgSpecTime /= count;
    avgMemUsage /= count;
    avgSpecies /= count;

    ss << "Average Times (seconds):\n";
    ss << "  Evaluation: " << avgEvalTime << "\n";
    ss << "  Mutation:   " << avgMutTime << "\n";
    ss << "  Crossover:  " << avgCrossTime << "\n";
    ss << "  Speciation: " << avgSpecTime << "\n\n";

    ss << "Memory and Population:\n";
    ss << "  Average Memory Usage (MB): " << avgMemUsage << "\n";
    ss << "  Average Active Species:    " << avgSpecies << "\n";
    ss << "  Total Innovations:         " << generationMetrics.back().innovations << "\n";
    ss << "  Total Evaluations:         " << generationMetrics.back().totalEvaluations << "\n";

    return ss.str();
}

void PerformanceMonitor::exportMetrics(const std::string& filename) {
    std::ofstream file(filename);
    if (!file) return;

    // Write CSV header
    file << "Generation,EvalTime,MutTime,CrossTime,SpecTime,TotalTime,"
         << "ActiveSpecies,TotalEvals,Innovations,MemoryUsage\n";

    // Write metrics
    for (size_t i = 0; i < generationMetrics.size(); ++i) {
        const auto& m = generationMetrics[i];
        file << i << ","
             << m.evaluationTime << ","
             << m.mutationTime << ","
             << m.crossoverTime << ","
             << m.speciationTime << ","
             << m.totalTime << ","
             << m.activeSpecies << ","
             << m.totalEvaluations << ","
             << m.innovations << ","
             << m.memoryUsage << "\n";
    }
}

}
}
