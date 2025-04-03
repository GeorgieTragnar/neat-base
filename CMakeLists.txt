cmake_minimum_required(VERSION 3.10)

# Setup vcpkg script with CMake (note: should be placed before project() call)
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/../vcpkg/scripts/buildsystems/vcpkg.cmake CACHE STRING "Vcpkg toolchain file")

project(NeatBase)

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

# Get all source files
file(GLOB_RECURSE SRC_FILES "${CMAKE_SOURCE_DIR}/*.cpp")
list(FILTER SRC_FILES EXCLUDE REGEX ".*tests/.*\\.cpp$")

# Define main.cpp path
set(MAIN_CPP "${CMAKE_SOURCE_DIR}/src/main.cpp")

# Remove main.cpp from library sources
list(REMOVE_ITEM SRC_FILES "${MAIN_CPP}")

# Create library target (all sources except main.cpp)
add_library(neat_lib STATIC ${SRC_FILES})
target_include_directories(neat_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(neat_lib PUBLIC 
    nlohmann_json::nlohmann_json
    spdlog::spdlog
)

# Create executable with main.cpp
add_executable(NeatBase "${MAIN_CPP}")
target_link_libraries(NeatBase PRIVATE neat_lib)

# Add test directory
add_subdirectory(tests)