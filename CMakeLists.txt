cmake_minimum_required(VERSION 3.10)

# Setup vcpkg script with CMake (note: should be placed before project() call)
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/../vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")

project(NeatBase)

# Add this near the top after project()
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "Debug flags")

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing
enable_testing()

# Find Google Test package
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}) 

# Find other dependencies
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

# Universal Library Code Generation Pipeline
include(version3/discover_all_types.cmake)
include(version3/generate_forward_declarations_from_types.cmake)
include(version3/generate_universal_operator_declarations.cmake)
include(version3/generate_template_friend_declarations.cmake)

# Define output files
set(DISCOVERED_TYPES_FILE 
    ${CMAKE_BINARY_DIR}/discovered_types.cmake)
set(GENERATED_DATA_FORWARDS_FILE 
    ${CMAKE_BINARY_DIR}/data_forward_declarations.inc)
set(GENERATED_OPERATOR_FORWARDS_FILE 
    ${CMAKE_BINARY_DIR}/operator_forward_declarations.inc)
set(GENERATED_OPERATOR_FRIENDS_FILE 
    ${CMAKE_BINARY_DIR}/operator_friend_declarations.inc)
set(GENERATED_TEMPLATE_FRIENDS_FILE 
    ${CMAKE_BINARY_DIR}/operator_template_friend_declarations.inc)

# Phase 1: Discover all types in the codebase
discover_all_types(
    ${CMAKE_SOURCE_DIR}/version3
    ${DISCOVERED_TYPES_FILE}
)

# Phase 2: Generate forward declarations from discovered types
generate_forward_declarations_from_types(
    ${DISCOVERED_TYPES_FILE}
    ${GENERATED_DATA_FORWARDS_FILE}
)

# Phase 3: Generate universal operator declarations
generate_universal_operator_declarations(
    ${CMAKE_SOURCE_DIR}/version3/operator
    ${DISCOVERED_TYPES_FILE}
    ${GENERATED_OPERATOR_FORWARDS_FILE}
    ${GENERATED_OPERATOR_FRIENDS_FILE}
)

# Phase 4: Generate template-specific friend declarations
generate_template_friend_declarations(
    ${CMAKE_SOURCE_DIR}/version3/operator
    ${DISCOVERED_TYPES_FILE}
    ${GENERATED_TEMPLATE_FRIENDS_FILE}
)

# Create custom commands to regenerate when source files change
add_custom_command(
    OUTPUT ${DISCOVERED_TYPES_FILE} ${GENERATED_DATA_FORWARDS_FILE} ${GENERATED_OPERATOR_FORWARDS_FILE} ${GENERATED_OPERATOR_FRIENDS_FILE} ${GENERATED_TEMPLATE_FRIENDS_FILE}
    COMMAND ${CMAKE_COMMAND} 
        -DCMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_SOURCE_DIR}/regenerate_all.cmake
    DEPENDS ${CMAKE_SOURCE_DIR}/version3/discover_all_types.cmake
            ${CMAKE_SOURCE_DIR}/version3/generate_forward_declarations_from_types.cmake
            ${CMAKE_SOURCE_DIR}/version3/generate_universal_operator_declarations.cmake
            ${CMAKE_SOURCE_DIR}/version3/generate_template_friend_declarations.cmake
    COMMENT "Regenerating all type discoveries and declarations"
    VERBATIM
)

# Get all source files
file(GLOB_RECURSE SRC_FILES "${CMAKE_SOURCE_DIR}/*.cpp")
list(FILTER SRC_FILES EXCLUDE REGEX ".*tests/.*\\.cpp$")

# Define main.cpp path
set(MAIN_CPP "${CMAKE_SOURCE_DIR}/version3/main.cpp")

# Remove main.cpp from library sources
list(REMOVE_ITEM SRC_FILES "${MAIN_CPP}")

# Create library target (all sources except main.cpp)
add_library(neat_lib STATIC ${SRC_FILES})
# Add the generated files as dependencies, not as sources
add_custom_target(generate_operator_declarations DEPENDS ${GENERATED_OPERATOR_FORWARDS_FILE} ${GENERATED_OPERATOR_FRIENDS_FILE})
add_dependencies(neat_lib generate_operator_declarations)

target_include_directories(neat_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}  # For version3 headers
    ${CMAKE_BINARY_DIR}  # For generated files
)
target_link_libraries(neat_lib PUBLIC 
    nlohmann_json::nlohmann_json
    spdlog::spdlog
)

# Create executable with main.cpp
add_executable(NeatBase "${MAIN_CPP}")
target_link_libraries(NeatBase PRIVATE neat_lib)

# Add test directory
# add_subdirectory(tests)