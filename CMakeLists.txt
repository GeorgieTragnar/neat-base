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

# Generate operator forward declarations
include(version3/generate_operator_forward_declarations.cmake)
include(version3/generate_operator_friend_declarations.cmake)

set(GENERATED_OPERATOR_FORWARDS_FILE 
    ${CMAKE_BINARY_DIR}/operator_forward_declarations.inc)
set(GENERATED_OPERATOR_FRIENDS_FILE 
    ${CMAKE_BINARY_DIR}/operator_friend_declarations.inc)

generate_operator_forward_declarations(
    ${CMAKE_SOURCE_DIR}/version3/operator
    ${GENERATED_OPERATOR_FORWARDS_FILE}
)

generate_operator_friend_declarations(
    ${CMAKE_SOURCE_DIR}/version3/operator
    ${GENERATED_OPERATOR_FRIENDS_FILE}
)

# Create a wrapper script for the generation
set(GENERATE_SCRIPT ${CMAKE_BINARY_DIR}/run_operator_generation.cmake)
file(WRITE ${GENERATE_SCRIPT}
"include(${CMAKE_SOURCE_DIR}/version3/generate_operator_forward_declarations.cmake)
include(${CMAKE_SOURCE_DIR}/version3/generate_operator_friend_declarations.cmake)
generate_operator_forward_declarations(
    \"${CMAKE_SOURCE_DIR}/version3/operator\"
    \"${GENERATED_OPERATOR_FORWARDS_FILE}\"
)
generate_operator_friend_declarations(
    \"${CMAKE_SOURCE_DIR}/version3/operator\"
    \"${GENERATED_OPERATOR_FRIENDS_FILE}\"
)")

# Create custom commands to regenerate when operator headers change
add_custom_command(
    OUTPUT ${GENERATED_OPERATOR_FORWARDS_FILE} ${GENERATED_OPERATOR_FRIENDS_FILE}
    COMMAND ${CMAKE_COMMAND} -P ${GENERATE_SCRIPT}
    DEPENDS ${CMAKE_SOURCE_DIR}/version3/generate_operator_forward_declarations.cmake
            ${CMAKE_SOURCE_DIR}/version3/generate_operator_friend_declarations.cmake
    COMMENT "Generating operator forward and friend declarations"
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