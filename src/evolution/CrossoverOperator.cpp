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
    
    
    // First inherit nodes from both parents to ensure they exist
    auto nodes = parent1.getNodes();
    nodes.insert(parent2.getNodes().begin(), parent2.getNodes().end());
    for (const auto& [id, type] : nodes) {
        child.addNode(id, type, false);  // Disable validation during construction
    }
    
    // Get genes from both parents - they should already be sorted by innovation number
    const auto& genes1 = parent1.getGenes();
    const auto& genes2 = parent2.getGenes();
    
    // Find max innovation number in each parent to identify excess genes
    int32_t maxInnov1 = genes1.empty() ? 0 : genes1.back().innovation;
    int32_t maxInnov2 = genes2.empty() ? 0 : genes2.back().innovation;
    int32_t smallerMaxInnov = std::min(maxInnov1, maxInnov2);
    
    LOG_DEBUG("Max innovation: Parent1 = {}, Parent2 = {}, Smaller = {}", 
              maxInnov1, maxInnov2, smallerMaxInnov);
    
    // Process genes in order of innovation number
    size_t i = 0, j = 0;
    while (i < genes1.size() && j < genes2.size()) {
        const core::Gene& gene1 = genes1[i];
        const core::Gene& gene2 = genes2[j];
        
        if (gene1.innovation == gene2.innovation) {
            // Matching gene - randomly inherit from either parent
            const core::Gene& selectedGene = 
                (inheritanceDist(rng) < config.matchingGeneInheritanceRate) 
                ? gene1 
                : gene2;
            
            core::Gene newGene = selectedGene;
            
            // Handle disabling: if gene is disabled in either parent, 75% chance to disable
            if (!gene1.enabled || !gene2.enabled) {
                newGene.enabled = !(inheritanceDist(rng) < 0.75);
            }
            
            child.addGene(newGene, false);
            i++;
            j++;
        }
        else if (gene1.innovation < gene2.innovation) {
            // Disjoint or excess gene from parent1
            // Inherit only if parent1 is fitter or equal fitness
            if (parent1IsFitter) {
                core::Gene newGene = gene1;
                child.addGene(newGene, false);
            }
            i++;
        }
        else {
            // Disjoint or excess gene from parent2
            // Inherit only if parent2 is fitter
            if (!parent1IsFitter) {
                core::Gene newGene = gene2;
                child.addGene(newGene, false);
            }
            j++;
        }
    }
    
    // Handle remaining genes from parent1 (must be excess genes)
    while (i < genes1.size()) {
        if (parent1IsFitter) {
            core::Gene newGene = genes1[i];
            child.addGene(newGene, false);
        }
        i++;
    }
    
    // Handle remaining genes from parent2 (must be excess genes)
    while (j < genes2.size()) {
        if (!parent1IsFitter) {
            core::Gene newGene = genes2[j];
            child.addGene(newGene, false);
        }
        j++;
    }
    
    LOG_INFO("Child created with {} nodes and {} genes", child.getNodes().size(), child.getGenes().size());
    
    child.rebuildNetwork();
    child.validate();
    
    return child;
}

}
}