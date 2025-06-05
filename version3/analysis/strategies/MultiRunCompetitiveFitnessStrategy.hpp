#pragma once
#include <functional>
#include <vector>
#include "../FitnessStrategy.hpp"

namespace Analysis {

enum class MultiChampionSelectionStrategy {
    BEST_CHAMPION,
    RANDOM_CHAMPION,
    ALL_CHAMPIONS
};

template<typename FitnessResultType>
struct MultiRunCompetitiveFitnessParams {
    std::function<FitnessResultType(const Phenotype&, const Phenotype&)> competitiveEvaluationFunction;
    int numFights;
    MultiChampionSelectionStrategy championSelection = MultiChampionSelectionStrategy::BEST_CHAMPION;
    std::function<FitnessResultType(const std::vector<FitnessResultType>&)> resolver;
};

template<typename FitnessResultType>
class MultiRunCompetitiveFitnessStrategy : public FitnessStrategy<FitnessResultType> {
public:
    explicit MultiRunCompetitiveFitnessStrategy(const MultiRunCompetitiveFitnessParams<FitnessResultType>& params)
        : _params(params) {}
    
    FitnessResultType evaluate(
        const Phenotype& candidate,
        const SpeciationControlUnit& speciation) const override {
        std::vector<FitnessResultType> results;
        results.reserve(_params.numFights);
        
        for (int i = 0; i < _params.numFights; ++i) {
            std::shared_ptr<const Phenotype> champion;
            
            switch (_params.championSelection) {
                case MultiChampionSelectionStrategy::BEST_CHAMPION:
                    champion = speciation.getBestChampion();
                    break;
                case MultiChampionSelectionStrategy::RANDOM_CHAMPION:
                    champion = speciation.getRandomChampion();
                    break;
                case MultiChampionSelectionStrategy::ALL_CHAMPIONS:
                    {
                        auto champions = speciation.getChampions();
                        champion = champions[i % champions.size()];
                    }
                    break;
            }
            
            results.push_back(_params.competitiveEvaluationFunction(candidate, *champion));
        }
        
        return _params.resolver(results);
    }

private:
    MultiRunCompetitiveFitnessParams<FitnessResultType> _params;
};

}