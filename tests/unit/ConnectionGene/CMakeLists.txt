# Include the generation functions
include(${CMAKE_CURRENT_SOURCE_DIR}/generate_connection_attributes_registration.cmake)

# Generate the field registration code
set(GENERATED_REGISTRATION_FILE 
    ${CMAKE_CURRENT_BINARY_DIR}/connection_gene_fields.inc)

generate_connection_attributes_registration(
    ${CMAKE_SOURCE_DIR}/version3/data/ConnectionGene.hpp
    ${GENERATED_REGISTRATION_FILE}
)

# Make the test file available to parent scope
set(CONNECTION_GENE_TEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ConnectionGene.cpp PARENT_SCOPE)

# Add include directory for generated file
set(CONNECTION_GENE_TEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

# Create custom targets to ensure regeneration when headers change
add_custom_command(
    OUTPUT ${GENERATED_REGISTRATION_FILE}
    COMMAND ${CMAKE_COMMAND} 
            -DSOURCE_FILE=${CMAKE_SOURCE_DIR}/version3/data/ConnectionGene.hpp
            -DOUTPUT_FILE=${GENERATED_REGISTRATION_FILE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_generation.cmake
    DEPENDS ${CMAKE_SOURCE_DIR}/version3/data/ConnectionGene.hpp
    COMMENT "Generating ConnectionGeneAttributes field registration"
)