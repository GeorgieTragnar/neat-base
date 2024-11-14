
#include "ReactComponents.hpp"

#include <sstream>

namespace neat {
namespace visualization {

static std::string generateNetworkViewer(const core::Genome& genome) {
	// Return a React component for network visualization
	// Using tailwind classes for styling
	std::ostringstream ss;
	ss << "export default function NetworkViewer() {\n";
	ss << "  return (\n";
	ss << "    <div className=\"w-full h-full flex flex-col\">\n";
	ss << "      <div className=\"flex-1 relative\">\n";
	ss << "        <svg className=\"w-full h-full\">\n";
	// Add SVG content here using NodeRenderer
	ss << "        </svg>\n";
	ss << "      </div>\n";
	ss << "      <div className=\"h-24 bg-gray-100 p-4\">\n";
	ss << "        <div className=\"flex justify-between\">\n";
	ss << "          <button className=\"px-4 py-2 bg-blue-500 text-white rounded\">\n";
	ss << "            Step\n";
	ss << "          </button>\n";
	ss << "          <button className=\"px-4 py-2 bg-green-500 text-white rounded\">\n";
	ss << "            Run\n";
	ss << "          </button>\n";
	ss << "          <button className=\"px-4 py-2 bg-red-500 text-white rounded\">\n";
	ss << "            Reset\n";
	ss << "          </button>\n";
	ss << "        </div>\n";
	ss << "      </div>\n";
	ss << "    </div>\n";
	ss << "  );\n";
	ss << "}\n";
	return ss.str();
}

static std::string generateStatsViewer(const StatsTracker& stats) {
	// Return a React component for statistics visualization
	std::ostringstream ss;
	ss << "export default function StatsViewer() {\n";
	ss << "  return (\n";
	ss << "    <div className=\"w-full h-full grid grid-cols-2 gap-4 p-4\">\n";
	ss << "      <div className=\"border rounded p-4\">\n";
	ss << "        <h3 className=\"text-lg font-semibold mb-2\">Fitness Over Time</h3>\n";
	// Add fitness plot component
	ss << "      </div>\n";
	ss << "      <div className=\"border rounded p-4\">\n";
	ss << "        <h3 className=\"text-lg font-semibold mb-2\">Species Distribution</h3>\n";
	// Add species plot component
	ss << "      </div>\n";
	ss << "    </div>\n";
	ss << "  );\n";
	ss << "}\n";
	return ss.str();
}

}
}
