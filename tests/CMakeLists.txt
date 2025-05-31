# Process subdirectories with their own CMakeLists.txt
set(ADDITIONAL_TEST_FILES "")
set(ADDITIONAL_INCLUDE_DIRS "")

# Check for unit test subdirectories
file(GLOB UNIT_SUBDIRS "${CMAKE_CURRENT_SOURCE_DIR}/unit/*")
foreach(SUBDIR ${UNIT_SUBDIRS})
    if(IS_DIRECTORY ${SUBDIR} AND EXISTS ${SUBDIR}/CMakeLists.txt)
        add_subdirectory(${SUBDIR})
        
        # Collect any test files from subdirectory
        if(DEFINED NODE_GENE_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${NODE_GENE_TEST_FILE})
        endif()
        if(DEFINED CONNECTION_GENE_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${CONNECTION_GENE_TEST_FILE})
        endif()
        if(DEFINED GENOME_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${GENOME_TEST_FILE})
        endif()
        
        # Collect any additional include directories
        if(DEFINED NODE_GENE_TEST_INCLUDE_DIR)
            list(APPEND ADDITIONAL_INCLUDE_DIRS ${NODE_GENE_TEST_INCLUDE_DIR})
        endif()
        if(DEFINED CONNECTION_GENE_TEST_INCLUDE_DIR)
            list(APPEND ADDITIONAL_INCLUDE_DIRS ${CONNECTION_GENE_TEST_INCLUDE_DIR})
        endif()
    endif()
endforeach()

# Check for integration test subdirectories
file(GLOB INTEGRATION_SUBDIRS "${CMAKE_CURRENT_SOURCE_DIR}/integration/*")
foreach(SUBDIR ${INTEGRATION_SUBDIRS})
    if(IS_DIRECTORY ${SUBDIR} AND EXISTS ${SUBDIR}/CMakeLists.txt)
        add_subdirectory(${SUBDIR})
        
        # Collect any test files from subdirectory
        if(DEFINED NODE_GENE_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${NODE_GENE_TEST_FILE})
        endif()
        if(DEFINED CONNECTION_GENE_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${CONNECTION_GENE_TEST_FILE})
        endif()
        if(DEFINED GENOME_TEST_FILE)
            list(APPEND ADDITIONAL_TEST_FILES ${GENOME_TEST_FILE})
        endif()
        
        # Collect any additional include directories
        if(DEFINED NODE_GENE_TEST_INCLUDE_DIR)
            list(APPEND ADDITIONAL_INCLUDE_DIRS ${NODE_GENE_TEST_INCLUDE_DIR})
        endif()
        if(DEFINED CONNECTION_GENE_TEST_INCLUDE_DIR)
            list(APPEND ADDITIONAL_INCLUDE_DIRS ${CONNECTION_GENE_TEST_INCLUDE_DIR})
        endif()
    endif()
endforeach()

# Test main file
file(GLOB_RECURSE UNIT_TEST_FILES "unit/*.cpp")
file(GLOB_RECURSE INTEGRATION_TEST_FILES "integration/*.cpp")

# Remove files that are handled by subdirectories
foreach(HANDLED_FILE ${ADDITIONAL_TEST_FILES})
    list(REMOVE_ITEM UNIT_TEST_FILES ${HANDLED_FILE})
    list(REMOVE_ITEM INTEGRATION_TEST_FILES ${HANDLED_FILE})
endforeach()

add_executable(neat_tests 
    main.cpp
    test_utilities.cpp
    ${UNIT_TEST_FILES}
    ${INTEGRATION_TEST_FILES}
    ${ADDITIONAL_TEST_FILES}
)

# Include directories
target_include_directories(neat_tests PUBLIC 
    ${CMAKE_SOURCE_DIR}
    ${ADDITIONAL_INCLUDE_DIRS}
)

# Add TESTING compile definition
target_compile_definitions(neat_tests PRIVATE TESTING)

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