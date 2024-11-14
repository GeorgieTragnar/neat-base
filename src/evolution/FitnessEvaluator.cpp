#include "FitnessEvaluator.hpp"
#include "core/Population.hpp"
#include "core/Genome.hpp"
#include "core/Species.hpp"
#include <algorithm>
#include <future>
#include <thread>

namespace neat {
namespace evolution {

FitnessEvaluator::FitnessEvaluator(const Config& config)
    : config(config)
{}

void FitnessEvaluator::evaluate(core::Population& population, const FitnessFunction& fitnessFunc) {
    if (config.useParallelEvaluation) {
        parallelEvaluate(population, fitnessFunc);
    } else {
        serialEvaluate(population, fitnessFunc);
    }
}

void FitnessEvaluator::parallelEvaluate(core::Population& population, const FitnessFunction& fitnessFunc) {
    const size_t numThreads = std::thread::hardware_concurrency();
    std::vector<std::future<void>> futures;
    
    for (size_t i = 0; i < population.size(); i += numThreads) {
        for (size_t j = 0; j < numThreads && i + j < population.size(); ++j) {
            futures.push_back(std::async(std::launch::async, [&, idx = i + j]() {
                auto& genome = population.getGenome(idx);
                double rawFitness = fitnessFunc(genome);
                genome.setFitness(processFitness(rawFitness));
            }));
        }
        
        for (auto& future : futures) {
            future.wait();
        }
        futures.clear();
    }
}

void FitnessEvaluator::serialEvaluate(core::Population& population, const FitnessFunction& fitnessFunc) {
    for (size_t i = 0; i < population.size(); ++i) {
        auto& genome = population.getGenome(i);
        double rawFitness = fitnessFunc(genome);
        genome.setFitness(processFitness(rawFitness));
    }
}

double FitnessEvaluator::processFitness(double rawFitness) const {
    if (!config.higherIsBetter) {
        rawFitness = config.targetFitness - rawFitness;
    }
    
    if (config.normalized) {
        rawFitness = std::max(0.0, rawFitness - config.fitnessThreshold);
        if (config.targetFitness > 0.0) {
            rawFitness /= config.targetFitness;
        }
    }
    
    return rawFitness;
}

void FitnessEvaluator::adjustFitness(core::Population& population) {
    // Adjust fitness based on species size to promote diversity
    for (size_t i = 0; i < population.size(); ++i) {
        auto& genome = population.getGenome(i);
        const auto& species = population.getSpecies();
        
        for (const auto& s : species) {
            if (s.getMembers().size() > 0) {
                double adjustedFitness = getAdjustedFitness(genome, s);
                genome.setAdjustedFitness(adjustedFitness);
            }
        }
    }
}

double FitnessEvaluator::getAdjustedFitness(const core::Genome& genome, const core::Species& species) const {
    // Explicit fitness sharing - divide by number of individuals in species
    return genome.getFitness() / species.getMembers().size();
}

}
}
