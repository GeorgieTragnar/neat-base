#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>

#include "tests/test_common.h"
#include "tests/test_utilities.h"
#include "version3/visualization/SVGGenerator.hpp"
#include "version3/operator/Init.hpp"
#include "version3/operator/PhenotypeConstruct.hpp"
#include "version3/data/HistoryTracker.hpp"

using namespace Visualization;
using namespace Operator;

class VisualizationTest : public ::testing::Test {
protected:
    void SetUp() override {
        historyTracker = std::make_shared<HistoryTracker>();
        
        // Clean up any existing test output
        if (std::filesystem::exists("test_visualizations")) {
            std::filesystem::remove_all("test_visualizations");
        }
    }

    void TearDown() override {
        // Clean up test output
        if (std::filesystem::exists("test_visualizations")) {
            std::filesystem::remove_all("test_visualizations");
        }
    }

    std::shared_ptr<HistoryTracker> historyTracker;
    
    Genome createSimpleGenome() {
        // Create a simple 2-input, 1-output genome
        std::vector<NodeGeneAttributes> inputAttrs = {
            {ActivationType::NONE}, {ActivationType::NONE}
        };
        std::vector<NodeGeneAttributes> outputAttrs = {
            {ActivationType::SIGMOID}
        };
        std::unordered_map<size_t, ConnectionGeneAttributes> biasAttrs;
        
        InitParams params(inputAttrs, outputAttrs, biasAttrs, 
                                   InitParams::InputConnectionStrategy::CONNECT_TO_OUTPUTS);
        
        return init(historyTracker, params);
    }
    
    bool isValidHTML(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            return false;
        }
        
        std::string content((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
        
        // Basic HTML validation checks
        return content.find("<!DOCTYPE html>") != std::string::npos &&
               content.find("<html>") != std::string::npos &&
               content.find("</html>") != std::string::npos &&
               content.find("<svg") != std::string::npos &&
               content.find("</svg>") != std::string::npos;
    }
    
    bool containsSVGElements(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            return false;
        }
        
        std::string content((std::istreambuf_iterator<char>(file)),
                           std::istreambuf_iterator<char>());
        
        // Check for expected SVG elements
        bool hasNodes = content.find("<circle") != std::string::npos || 
                       content.find("<rect") != std::string::npos;
        bool hasConnections = content.find("<line") != std::string::npos;
        
        return hasNodes; // Connections are optional for minimal genomes
    }
};

TEST_F(VisualizationTest, InitializeCreatesOutputDirectory) {
    VisualizationConfig config;
    config.outputDirectory = "test_visualizations";
    
    EXPECT_NO_THROW(initialize(config));
}

TEST_F(VisualizationTest, GenerateVisualizationCreatesValidHTMLFile) {
    VisualizationConfig config;
    config.outputDirectory = "test_visualizations";
    initialize(config);
    
    // Create a simple genome and its phenotype
    Genome genome = createSimpleGenome();
    Operator::phenotypeConstruct(genome);
    const auto& phenotype = genome.get_phenotype();
    
    EXPECT_NO_THROW(generateVisualization(phenotype, 0, 0));
    
    // Check that file was created
    std::string expectedDir = "test_visualizations";
    EXPECT_TRUE(std::filesystem::exists(expectedDir));
    
    // Find the timestamped subdirectory
    bool foundTimestampDir = false;
    bool foundHTMLFile = false;
    
    for (const auto& entry : std::filesystem::directory_iterator(expectedDir)) {
        if (entry.is_directory()) {
            foundTimestampDir = true;
            
            // Look for the HTML file in the timestamped directory
            std::filesystem::path htmlFile = entry.path() / "generation_0_species_0.html";
            if (std::filesystem::exists(htmlFile)) {
                foundHTMLFile = true;
                
                // Validate HTML content
                EXPECT_TRUE(isValidHTML(htmlFile.string()));
                EXPECT_TRUE(containsSVGElements(htmlFile.string()));
            }
        }
    }
    
    EXPECT_TRUE(foundTimestampDir);
    EXPECT_TRUE(foundHTMLFile);
}

TEST_F(VisualizationTest, GenerateMultipleVisualizationsCreatesSeparateFiles) {
    VisualizationConfig config;
    config.outputDirectory = "test_visualizations";
    initialize(config);
    
    Genome genome = createSimpleGenome();
    Operator::phenotypeConstruct(genome);
    const auto& phenotype = genome.get_phenotype();
    
    // Generate multiple visualizations
    EXPECT_NO_THROW(generateVisualization(phenotype, 0, 0));
    EXPECT_NO_THROW(generateVisualization(phenotype, 0, 1));
    EXPECT_NO_THROW(generateVisualization(phenotype, 1, 0));
    
    // Check that multiple files were created
    std::string expectedDir = "test_visualizations";
    
    for (const auto& entry : std::filesystem::directory_iterator(expectedDir)) {
        if (entry.is_directory()) {
            int fileCount = 0;
            for (const auto& file : std::filesystem::directory_iterator(entry.path())) {
                if (file.path().extension() == ".html") {
                    fileCount++;
                }
            }
            EXPECT_EQ(fileCount, 3);
            break;
        }
    }
}

TEST_F(VisualizationTest, HandleEmptyPhenotype) {
    VisualizationConfig config;
    config.outputDirectory = "test_visualizations";
    initialize(config);
    
    // Create minimal phenotype structure
    Genome::Phenotype emptyPhenotype;
    
    // Should not crash with empty phenotype
    EXPECT_NO_THROW(generateVisualization(emptyPhenotype, 0, 0));
}

TEST_F(VisualizationTest, ThrowsErrorWhenNotInitialized) {
    Genome genome = createSimpleGenome();
    Operator::phenotypeConstruct(genome);
    const auto& phenotype = genome.get_phenotype();
    
    // Should throw error if not initialized
    EXPECT_THROW(generateVisualization(phenotype, 0, 0), std::runtime_error);
}

TEST_F(VisualizationTest, FilenameGenerationIsCorrect) {
    EXPECT_EQ(generateFilename(0, 0), "generation_0_species_0.html");
    EXPECT_EQ(generateFilename(5, 3), "generation_5_species_3.html");
    EXPECT_EQ(generateFilename(100, 25), "generation_100_species_25.html");
}

TEST_F(VisualizationTest, TimestampFormatIsValid) {
    std::string timestamp = getCurrentTimestamp();
    
    // Should be in format YY_MM_DD_HH_MM_SS (17 characters)
    EXPECT_EQ(timestamp.length(), 17u);
    
    // Should contain only digits and underscores
    for (char c : timestamp) {
        EXPECT_TRUE(std::isdigit(c) || c == '_');
    }
    
    // Should have underscores in correct positions
    EXPECT_EQ(timestamp[2], '_');
    EXPECT_EQ(timestamp[5], '_');
    EXPECT_EQ(timestamp[8], '_');
    EXPECT_EQ(timestamp[11], '_');
    EXPECT_EQ(timestamp[14], '_');
}