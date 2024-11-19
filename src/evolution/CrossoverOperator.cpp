#include "CrossoverOperator.hpp"
#include "core/Gene.hpp"
#include "core/Genome.hpp"
#include "core/InnovationTracker.hpp"

#include "logger/Logger.hpp"
static auto logger = LOGGER("evolution::CrossoverOperator");

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
    
    LOG_DEBUG("Starting crossover");
    LOG_DEBUG("Parent 1 fitness: {}, nodes: {}, genes: {}", parent1.getFitness(), parent1.getNodes().size(), parent1.getGenes().size());
    LOG_DEBUG("Parent 2 fitness: {}, nodes: {}, genes: {}", parent2.getFitness(), parent2.getNodes().size(), parent2.getGenes().size());
        
    // Create child genome with same config as parent
    core::Genome child(parent1.getConfig());
    child.getGenes().clear();
    
    // Get fitness values to determine which parent is more fit
    const double fitness1 = parent1.getFitness();
    const double fitness2 = parent2.getFitness();
    const bool parent1IsFitter = fitness1 >= fitness2;
    
    LOG_DEBUG("Parent 1 {} fitter", (parent1IsFitter ? "is" : "is not"));
    
    // Create inheritance map of matching and disjoint genes
    auto inheritanceMap = createInheritanceMap(parent1, parent2);
    
    // First inherit nodes from both parents to ensure they exist
    auto nodes = parent1.getNodes();
    nodes.insert(parent2.getNodes().begin(), parent2.getNodes().end());
    for (const auto& [id, type] : nodes) {
        child.addNode(id, type, false);  // Disable validation during construction
    }
    
    // Inherit genes according to rules
    LOG_DEBUG("Processing {} genes", inheritanceMap.size());
    
    
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
            
            child.addGene(newGene, false);  // Disable validation during construction
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
            
            child.addGene(newGene, false);  // Disable validation during construction
        }
    }
    
    LOG_INFO("Child created with {} nodes and {} genes", child.getNodes().size(), child.getGenes().size());
    
    child.rebuildNetwork();
    child.validate();
    
    return child;
}

}
}