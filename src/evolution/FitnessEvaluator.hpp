// fitness.hpp
#pragma once

#include <functional>

namespace neat {
namespace core {
class Genome;
class Population;
class Species;
}
namespace evolution {

class FitnessEvaluator {
public:
    struct Config {
        bool normalized = true;
        double targetFitness = 1.0;
        bool higherIsBetter = true;
        double fitnessThreshold = 0.0;
        bool useParallelEvaluation = true;
    };

    using FitnessFunction = std::function<double(const core::Genome&)>;
    
    explicit FitnessEvaluator(const Config& config);
    
    void evaluate(core::Population& population, const FitnessFunction& fitnessFunc);
    void adjustFitness(core::Population& population);
    double getAdjustedFitness(const core::Genome& genome, const core::Species& species) const;

private:
    double processFitness(double rawFitness) const;
    void parallelEvaluate(core::Population& population, const FitnessFunction& fitnessFunc);
    void serialEvaluate(core::Population& population, const FitnessFunction& fitnessFunc);

    Config config;
};

}
}
