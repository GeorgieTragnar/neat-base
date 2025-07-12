# PopulationIntegration Integration Test CMakeLists.txt
# This test will be included in the main neat_tests executable via the parent test system

# Set the test file for inclusion in the main test executable
set(POPULATION_INTEGRATION_TEST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/CoreGenomePropagation.cpp" PARENT_SCOPE)