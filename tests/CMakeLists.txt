# Test main file
file(GLOB_RECURSE UNIT_TEST_FILES "unit/*.cpp")
file(GLOB_RECURSE INTEGRATION_TEST_FILES "integration/*.cpp")

add_executable(neat_tests 
    main.cpp
    ${UNIT_TEST_FILES}
    ${INTEGRATION_TEST_FILES}
)

# Include directories
target_include_directories(neat_tests PUBLIC 
    ${CMAKE_SOURCE_DIR}
)

# Link with your NEAT library and GTest
target_link_libraries(neat_tests
    neat_lib
    GTest::GTest
    GTest::Main
    pthread
)

# Discover tests
include(GoogleTest)
gtest_discover_tests(neat_tests)