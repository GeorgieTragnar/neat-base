# Include the generation functions
include(${CMAKE_CURRENT_SOURCE_DIR}/generate_node_attributes_registration.cmake)

# Generate the field registration code
set(GENERATED_REGISTRATION_FILE 
    ${CMAKE_CURRENT_BINARY_DIR}/node_gene_fields.inc)

generate_node_attributes_registration(
    ${CMAKE_SOURCE_DIR}/version3/data/NodeGene.hpp
    ${GENERATED_REGISTRATION_FILE}
)

# Make the test file available to parent scope
set(NODE_GENE_TEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/NodeGene.cpp PARENT_SCOPE)

# Add include directory for generated file
set(NODE_GENE_TEST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

# Create custom targets to ensure regeneration when headers change
add_custom_command(
    OUTPUT ${GENERATED_REGISTRATION_FILE}
    COMMAND ${CMAKE_COMMAND} 
            -DSOURCE_FILE=${CMAKE_SOURCE_DIR}/version3/data/NodeGene.hpp
            -DOUTPUT_FILE=${GENERATED_REGISTRATION_FILE}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_generation.cmake
    DEPENDS ${CMAKE_SOURCE_DIR}/version3/data/NodeGene.hpp
    COMMENT "Generating NodeGeneAttributes field registration"
)