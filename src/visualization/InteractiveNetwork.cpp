#include "InteractiveNetwork.hpp"
#include "core/Genome.hpp"
#include "network/Network.hpp"
#include <sstream>
#include <chrono>
#include <thread>
#include <nlohmann/json.hpp>

namespace neat {
namespace visualization {

InteractiveNetwork::InteractiveNetwork(const core::Genome& genome, const Config& config)
    : genome(genome)
    , config(config)
{
    reset();
}

void InteractiveNetwork::start() {
    network::Network::Config netConfig;
    network::Network net(netConfig);
    
    for (const auto& [id, type] : genome.getNodes()) {
        net.addNode(id, type);
    }
    
    for (const auto& gene : genome.getGenes()) {
        if (gene.enabled) {
            net.addConnection(gene.inputNode, gene.outputNode, gene.weight);
        }
    }
    
    while (true) {
        step();
        if (config.updateInterval > 0) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(config.updateInterval)
            );
        }
    }
}

void InteractiveNetwork::updateInputs(const std::vector<double>& inputs) {
    int inputIdx = 0;
    for (const auto& [id, type] : genome.getNodes()) {
        if (type == core::ENodeType::INPUT && inputIdx < inputs.size()) {
            nodeValues[id] = inputs[inputIdx++];
        }
    }
    
    std::vector<double> currentValues;
    for (const auto& [id, value] : nodeValues) {
        currentValues.push_back(value);
    }
    signalHistory.push_back(currentValues);
}

void InteractiveNetwork::step() {
    network::Network::Config netConfig;
    network::Network net(netConfig);
    
    for (const auto& [id, type] : genome.getNodes()) {
        net.addNode(id, type);
    }
    
    for (const auto& gene : genome.getGenes()) {
        if (gene.enabled) {
            net.addConnection(gene.inputNode, gene.outputNode, gene.weight);
        }
    }
    
    std::vector<double> inputs;
    for (const auto& [id, type] : genome.getNodes()) {
        if (type == core::ENodeType::INPUT) {
            inputs.push_back(nodeValues[id]);
        }
    }
    
    auto outputs = net.activate(inputs);
    
    int outputIdx = 0;
    for (const auto& [id, type] : genome.getNodes()) {
        if (type == core::ENodeType::OUTPUT && outputIdx < outputs.size()) {
            nodeValues[id] = outputs[outputIdx++];
        }
    }
    
    std::vector<double> currentValues;
    for (const auto& [id, value] : nodeValues) {
        currentValues.push_back(value);
    }
    signalHistory.push_back(currentValues);
}

void InteractiveNetwork::reset() {
    nodeValues.clear();
    signalHistory.clear();
    
    for (const auto& [id, type] : genome.getNodes()) {
        nodeValues[id] = 0.0;
    }
}

void InteractiveNetwork::onNodeClick(int32_t nodeId) {
    if (nodeValues.find(nodeId) != nodeValues.end()) {
        nodeValues[nodeId] = nodeValues[nodeId] > 0.5 ? 0.0 : 1.0;
    }
}

void InteractiveNetwork::onConnectionClick(int32_t fromId, int32_t toId) {
    // Handle connection click event
    // Note: Genome is const, so we can only observe connections
}

std::string InteractiveNetwork::generateReactComponent() {
    nlohmann::json nodeValuesJson;
    for (const auto& [id, value] : nodeValues) {
        nodeValuesJson[std::to_string(id)] = value;
    }
    
    std::ostringstream ss;
    ss << "export default function NetworkViewer({ onNodeClick, onConnectionClick }) {\n"
       << "  const [values, setValues] = useState(" << nodeValuesJson.dump() << ");\n"
       << "  const [history, setHistory] = useState([]);\n\n"
       << "  useEffect(() => {\n"
       << "    const interval = setInterval(() => {\n"
       << "      step();\n"
       << "    }, " << config.updateInterval << ");\n"
       << "    return () => clearInterval(interval);\n"
       << "  }, []);\n\n";
    
    // Add network visualization JSX
    ss << "  return (\n"
       << "    <div className=\"network-viewer\">\n"
       << "      <div className=\"nodes\">\n";
    
    for (const auto& [id, type] : genome.getNodes()) {
        ss << "        <Node id=\"" << id << "\" type=\"" << static_cast<int>(type) 
           << "\" value={values[" << id << "]} onClick={onNodeClick} />\n";
    }
    
    ss << "      </div>\n"
       << "      <div className=\"connections\">\n";
       
    for (const auto& gene : genome.getGenes()) {
        if (gene.enabled) {
            ss << "        <Connection from=\"" << gene.inputNode 
               << "\" to=\"" << gene.outputNode 
               << "\" weight=\"" << gene.weight 
               << "\" onClick={onConnectionClick} />\n";
        }
    }
    
    ss << "      </div>\n"
       << "    </div>\n"
       << "  );\n"
       << "}\n";
    
    return ss.str();
}

}
}
