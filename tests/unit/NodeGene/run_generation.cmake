# Script to run the generation function
# This is needed because custom_command can't directly call functions

include(${CMAKE_CURRENT_LIST_DIR}/generate_node_attributes_registration.cmake)

generate_node_attributes_registration(${SOURCE_FILE} ${OUTPUT_FILE})