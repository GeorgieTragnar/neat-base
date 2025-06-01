
#include "Init.hpp"

using namespace Operator;

InitializationParams::InitializationParams(const std::vector<NodeGeneAttributes>& inputAttributes, 
	const std::vector<NodeGeneAttributes>& outputAttributes,
	const std::unordered_map<size_t, ConnectionGeneAttributes>& biasAttributes,
	const InputConnectionStrategy& inputStrategy)
	: _inputAttributes(inputAttributes)
	, _outputAttributes(outputAttributes)
	, _biasAttributes(biasAttributes)
	, _inputStrategy(inputStrategy)
{

}

Genome Operator::init(std::unique_ptr<HistoryTracker>& historyTracker, const InitializationParams& params)
{
    const auto& inputAttributes = params._inputAttributes;
    const auto& outputAttributes = params._outputAttributes;
    const auto& biasAttributes = params._biasAttributes;
    const auto inputStrategy = params._inputStrategy;

    GenomeParams genomeParams;

    const size_t numInputs = inputAttributes.size();
    for (size_t i = 0; i < numInputs; i++) {
        uint32_t nodeID = historyTracker->get_input(i);
        genomeParams._nodeHistoryIDs.push_back(nodeID);
        genomeParams._nodeTypes.push_back(NodeType::INPUT);
        genomeParams._nodeAttributes.push_back(inputAttributes[i]);
    }
    
    // Add bias node gene
    const uint32_t biasNodeID = historyTracker->get_bias();
    genomeParams._nodeHistoryIDs.push_back(biasNodeID);
    genomeParams._nodeTypes.push_back(NodeType::BIAS);
    genomeParams._nodeAttributes.push_back(NodeGeneAttributes{ActivationType::NONE}); // Bias has no activation
    
    // Add node genes for outputs
    const size_t numOutputs = outputAttributes.size();
    for (size_t i = 0; i < numOutputs; i++) {
        uint32_t nodeID = historyTracker->get_output(i);
        genomeParams._nodeHistoryIDs.push_back(nodeID);
        genomeParams._nodeTypes.push_back(NodeType::OUTPUT);
        genomeParams._nodeAttributes.push_back(outputAttributes[i]);
    }
    
    // Create connections according to the strategy
    if (inputStrategy == InitializationParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS) {
        // Connect all inputs to all outputs
        for (size_t inputIdx = 0; inputIdx < numInputs; inputIdx++) {
            uint32_t inputNodeID = historyTracker->get_input(inputIdx);
            
            for (size_t outputIdx = 0; outputIdx < numOutputs; outputIdx++) {
                uint32_t outputNodeID = historyTracker->get_output(outputIdx);
                
                // Get innovation number for this connection
                uint32_t connectionID = historyTracker->get_connection(inputNodeID, outputNodeID);
                
                // Create connection
                genomeParams._connectionHistoryIDs.push_back(connectionID);
                genomeParams._sourceNodeHistoryIDs.push_back(inputNodeID);
                genomeParams._targetNodeHistoryIDs.push_back(outputNodeID);
                
                // Default connection attributes
                ConnectionGeneAttributes connAttr;
                connAttr.weight = (static_cast<float>(rand()) / RAND_MAX) * 4.0f - 2.0f; // Random weight between -2 and 2
                connAttr.enabled = true;
                
                genomeParams._connectionAttributes.push_back(connAttr);
            }
        }
        
        // Connect bias to outputs with provided attributes
        uint32_t biasID = historyTracker->get_bias();
        
        for (size_t outputIdx = 0; outputIdx < numOutputs; outputIdx++) {
            uint32_t outputNodeID = historyTracker->get_output(outputIdx);
            
            // Get innovation number for this connection
            uint32_t connectionID = historyTracker->get_connection(biasID, outputNodeID);
            
            // Create connection
            genomeParams._connectionHistoryIDs.push_back(connectionID);
            genomeParams._sourceNodeHistoryIDs.push_back(biasID);
            genomeParams._targetNodeHistoryIDs.push_back(outputNodeID);
            
            // Use bias connection attributes if provided, otherwise default
            ConnectionGeneAttributes connAttr;
            auto biasAttrIt = biasAttributes.find(outputIdx);
            if (biasAttrIt != biasAttributes.end()) {
                connAttr = biasAttrIt->second;
            } else {
                connAttr.weight = 1.0f;
                connAttr.enabled = true;
            }
            
            genomeParams._connectionAttributes.push_back(connAttr);
        }
    }
    
    // Create and return the genome with the constructed parameters
    return Genome(genomeParams);
}
