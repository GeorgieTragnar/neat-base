#include "NetworkExecution.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "../logger/Logger.hpp"

namespace Operator {

NetworkExecutionParams::NetworkExecutionParams(bool debugOutput)
    : _debugOutput(debugOutput) {
}

std::vector<double> networkExecution(
    const Genome::Phenotype& phenotype,
    const std::vector<double>& inputs,
    const NetworkExecutionParams& params
) {
    static auto logger = LOGGER("operator.NetworkExecution");
    
    const auto& connections = phenotype._orderedConnections;
    const auto& nodes = phenotype._nodeGeneAttributes;
    const auto& inputIndices = phenotype._inputIndices;
    const auto& outputIndices = phenotype._outputIndices;
    const size_t biasIndex = phenotype._biasIndex;
    
    // Size validation
    assert(inputs.size() == inputIndices.size() && "Input size must match number of input nodes");
    assert(!outputIndices.empty() && "Network must have at least one output node");
    
    if (params._debugOutput) {
        LOG_DEBUG("NetworkExecution: {} nodes, {} connections, {} inputs, {} outputs", 
                 nodes.size(), connections.size(), inputIndices.size(), outputIndices.size());
    }
    
    // Initialization Phase: Clear working memory and set initial values
    std::vector<double> nodeValues(nodes.size(), 0.0);
    
    // Set input values
    for (size_t i = 0; i < inputs.size(); ++i) {
        nodeValues[inputIndices[i]] = inputs[i];
        if (params._debugOutput) {
            LOG_DEBUG("Set input node[{}] = {:.3f}", inputIndices[i], inputs[i]);
        }
    }
    
    // Set bias value if bias node exists
    if (biasIndex != SIZE_MAX && biasIndex < nodeValues.size()) {
        nodeValues[biasIndex] = 1.0;
        if (params._debugOutput) {
            LOG_DEBUG("Set bias node[{}] = 1.0", biasIndex);
        }
    }
    
    // Dependency Analysis Phase: Build dependency graph and compute evaluation order
    std::vector<std::vector<size_t>> dependents(nodes.size()); // dependents[i] = nodes that depend on node i
    std::vector<size_t> incomingCount(nodes.size(), 0);       // number of incoming dependencies
    
    // Group connections by target node for efficient processing
    std::unordered_map<size_t, std::vector<size_t>> connectionsByTarget;
    
    for (size_t connIdx = 0; connIdx < connections.size(); ++connIdx) {
        const auto& conn = connections[connIdx];
        if (conn._connectionGeneAttribute.enabled) {
            size_t sourceIdx = conn._sourceNodeIndex;
            size_t targetIdx = conn._targetNodeIndex;
            
            // Build dependency graph
            dependents[sourceIdx].push_back(targetIdx);
            incomingCount[targetIdx]++;
            
            // Group connections by target
            connectionsByTarget[targetIdx].push_back(connIdx);
        }
    }
    
    // Compute topological order using Kahn's algorithm
    std::queue<size_t> readyQueue;
    std::vector<size_t> evaluationOrder;
    
    // Start with nodes that have no dependencies (inputs and bias)
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (incomingCount[i] == 0) {
            readyQueue.push(i);
        }
    }
    
    while (!readyQueue.empty()) {
        size_t currentNode = readyQueue.front();
        readyQueue.pop();
        evaluationOrder.push_back(currentNode);
        
        // Update dependents
        for (size_t dependent : dependents[currentNode]) {
            incomingCount[dependent]--;
            if (incomingCount[dependent] == 0) {
                readyQueue.push(dependent);
            }
        }
    }
    
    // Validate that all nodes were included (no cycles)
    assert(evaluationOrder.size() == nodes.size() && "Topological sort failed - network has cycles");
    
    if (params._debugOutput) {
        LOG_DEBUG("Evaluation order computed: {} nodes", evaluationOrder.size());
    }
    
    // Create sets for efficient pre-set node identification
    std::unordered_set<size_t> inputSet(inputIndices.begin(), inputIndices.end());
    bool hasBias = (biasIndex != SIZE_MAX && biasIndex < nodeValues.size());
    
    // Network Execution Phase: Process nodes in evaluation order
    for (size_t nodeIdx : evaluationOrder) {
        // Skip pre-set nodes (inputs and bias)
        if (inputSet.find(nodeIdx) != inputSet.end() || (hasBias && nodeIdx == biasIndex)) {
            if (params._debugOutput) {
                LOG_DEBUG("Skipping pre-set node[{}] = {:.3f}", nodeIdx, nodeValues[nodeIdx]);
            }
            continue;
        }
        
        // Accumulate weighted inputs from all incoming connections
        double accumulatedInput = 0.0;
        
        if (connectionsByTarget.find(nodeIdx) != connectionsByTarget.end()) {
            for (size_t connIdx : connectionsByTarget[nodeIdx]) {
                const auto& conn = connections[connIdx];
                if (conn._connectionGeneAttribute.enabled) {
                    double contribution = nodeValues[conn._sourceNodeIndex] * conn._connectionGeneAttribute.weight;
                    accumulatedInput += contribution;
                    
                    if (params._debugOutput) {
                        LOG_DEBUG("Connection: node[{}]({:.3f}) -> node[{}], weight={:.3f}, contribution={:.3f}",
                                 conn._sourceNodeIndex, nodeValues[conn._sourceNodeIndex], 
                                 nodeIdx, conn._connectionGeneAttribute.weight, contribution);
                    }
                }
            }
        }
        
        // Apply activation function
        double activatedValue;
        switch (nodes[nodeIdx].activationType) {
            case ActivationType::NONE:
                activatedValue = accumulatedInput;
                break;
                
            case ActivationType::SIGMOID: {
                // Clamp input to prevent overflow/underflow
                double clampedInput = std::max(-50.0, std::min(50.0, accumulatedInput));
                activatedValue = 1.0 / (1.0 + std::exp(-clampedInput));
                break;
            }
            
            case ActivationType::TANH: {
                // Clamp input to prevent overflow/underflow  
                double clampedInput = std::max(-50.0, std::min(50.0, accumulatedInput));
                activatedValue = std::tanh(clampedInput);
                break;
            }
            
            case ActivationType::RELU:
                activatedValue = std::max(0.0, accumulatedInput);
                break;
                
            case ActivationType::LEAKY_RELU: {
                const double alpha = 0.01; // Standard leaky coefficient
                activatedValue = accumulatedInput > 0.0 ? accumulatedInput : alpha * accumulatedInput;
                break;
            }
            
            case ActivationType::STEP:
                activatedValue = accumulatedInput > 0.0 ? 1.0 : 0.0;
                break;
                
            default:
                // Fallback to linear for unknown activation types
                activatedValue = accumulatedInput;
                break;
        }

        double oldValue = nodeValues[nodeIdx];
        nodeValues[nodeIdx] = activatedValue;
        
        if (params._debugOutput) {
            LOG_DEBUG("Activate node[{}]: input={:.3f} -> output={:.3f} (was {:.3f})", 
                     nodeIdx, accumulatedInput, nodeValues[nodeIdx], oldValue);
        }
    }
    
    // Output Extraction Phase: Copy output values
    std::vector<double> outputs;
    outputs.reserve(outputIndices.size());
    
    for (size_t outputIdx : outputIndices) {
        outputs.push_back(nodeValues[outputIdx]);
        if (params._debugOutput) {
            LOG_DEBUG("Output node[{}] = {:.3f}", outputIdx, nodeValues[outputIdx]);
        }
    }
    
    if (params._debugOutput) {
        LOG_DEBUG("NetworkExecution completed: {} outputs generated", outputs.size());
    }
    
    return outputs;
}

} // namespace Operator