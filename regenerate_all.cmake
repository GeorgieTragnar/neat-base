# Regeneration script for all type discoveries and declarations
# This script is called by CMake custom commands when source files change

# Include all the generation functions
include(${CMAKE_SOURCE_DIR}/version3/discover_all_types.cmake)
include(${CMAKE_SOURCE_DIR}/version3/generate_forward_declarations_from_types.cmake)
include(${CMAKE_SOURCE_DIR}/version3/generate_universal_operator_declarations.cmake)
include(${CMAKE_SOURCE_DIR}/version3/generate_template_friend_declarations.cmake)

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

message(STATUS "Regenerating all type discoveries and declarations...")

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
    ${CMAKE_SOURCE_DIR}/version3
    ${DISCOVERED_TYPES_FILE}
    ${GENERATED_OPERATOR_FORWARDS_FILE}
    ${GENERATED_OPERATOR_FRIENDS_FILE}
)

# Phase 4: Generate template-specific friend declarations
generate_template_friend_declarations(
    ${CMAKE_SOURCE_DIR}/version3
    ${DISCOVERED_TYPES_FILE}
    ${GENERATED_TEMPLATE_FRIENDS_FILE}
)

message(STATUS "All type discoveries and declarations regenerated successfully!")