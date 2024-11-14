#include "CrossoverOperator.hpp"
#include "core/Gene.hpp"
#include "core/Genome.hpp"
#include "core/InnovationTracker.hpp"

namespace neat {
namespace evolution {

CrossoverOperator::CrossoverOperator(const Config& config)
    : config(config)
    , rng(std::random_device{}())
    , inheritanceDist(0.0, 1.0)
{}

CrossoverOperator::InheritanceMap CrossoverOperator::createInheritanceMap(
    const core::Genome& parent1, 
    const core::Genome& parent2) {
    
    InheritanceMap inheritanceMap;
    const auto& genes1 = parent1.getGenes();
    const auto& genes2 = parent2.getGenes();
    
    // Create mapping of innovation numbers to genes from both parents
    for (const auto& gene : genes1) {
        inheritanceMap[gene.innovation].first = std::make_shared<core::Gene>(gene);
    }
    
    for (const auto& gene : genes2) {
        inheritanceMap[gene.innovation].second = std::make_shared<core::Gene>(gene);
    }
    
    return inheritanceMap;
}

core::Genome CrossoverOperator::crossover(
    const core::Genome& parent1,
    const core::Genome& parent2) {
    
    // Create child genome with same config as parent
    core::Genome child(parent1.getConfig());
    child.getGenes().clear();
    
    // Get fitness values to determine which parent is more fit
    const double fitness1 = parent1.getFitness();
    const double fitness2 = parent2.getFitness();
    const bool parent1IsFitter = fitness1 >= fitness2;
    
    // Create inheritance map of matching and disjoint genes
    auto inheritanceMap = createInheritanceMap(parent1, parent2);
    
    // Inherit genes according to rules
    for (const auto& [innovation, genePair] : inheritanceMap) {
        const auto& [gene1, gene2] = genePair;
        
        // Both parents have this gene (matching)
        if (gene1 && gene2) {
            // Randomly inherit from either parent
            const auto& selectedGene = 
                (inheritanceDist(rng) < config.matchingGeneInheritanceRate) 
                ? gene1 
                : gene2;
                
            auto newGene = *selectedGene;
            
            // Check if gene should be re-enabled
            if (!newGene.enabled && 
                config.inheritDisabledGenes &&
                inheritanceDist(rng) < config.disabledGeneReenableRate) {
                newGene.enabled = true;
            }
            
            child.addGene(newGene);
        }
        // Disjoint/excess gene - inherit from fitter parent only
        else if ((gene1 && parent1IsFitter) || (gene2 && !parent1IsFitter)) {
            const auto& selectedGene = gene1 ? gene1 : gene2;
            auto newGene = *selectedGene;
            
            if (!newGene.enabled && 
                config.inheritDisabledGenes && 
                inheritanceDist(rng) < config.disabledGeneReenableRate) {
                newGene.enabled = true;
            }
            
            child.addGene(newGene);
        }
    }
    
    // Inherit nodes from both parents
    auto nodes = parent1.getNodes();
    nodes.insert(parent2.getNodes().begin(), parent2.getNodes().end());
    for (const auto& [id, type] : nodes) {
        child.addNode(id, type);
    }
    
    return child;
}

}
}