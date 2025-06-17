#pragma once
#include <functional>
#include "../FitnessStrategy.hpp"

namespace Analysis {

enum class ChampionSelectionStrategy {
    BEST_CHAMPION,
    RANDOM_CHAMPION
};

template<typename FitnessResultType>
struct SingleCompetitiveFitnessParams {
    std::function<FitnessResultType(const Phenotype&, const Phenotype&)> competitiveEvaluationFunction;
    ChampionSelectionStrategy championSelection = ChampionSelectionStrategy::BEST_CHAMPION;
};

template<typename FitnessResultType>
class SingleCompetitiveFitnessStrategy : public FitnessStrategy<FitnessResultType> {
public:
    explicit SingleCompetitiveFitnessStrategy(const SingleCompetitiveFitnessParams<FitnessResultType>& params)
        : _params(params) {}
    
    FitnessResultType evaluate(
        const Phenotype& candidate,
        const SpeciationControlUnit& speciation) const override {
        std::shared_ptr<const Phenotype> champion;
        
        switch (_params.championSelection) {
            case ChampionSelectionStrategy::BEST_CHAMPION:
                champion = speciation.getBestChampion();
                break;
            case ChampionSelectionStrategy::RANDOM_CHAMPION:
                champion = speciation.getRandomChampion();
                break;
        }
        
        return _params.competitiveEvaluationFunction(candidate, *champion);
    }

private:
    SingleCompetitiveFitnessParams<FitnessResultType> _params;
};

}