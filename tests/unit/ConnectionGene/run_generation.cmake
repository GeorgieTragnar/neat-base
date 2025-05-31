# This script is run by the custom command to generate field registration
include(${CMAKE_CURRENT_LIST_DIR}/generate_connection_attributes_registration.cmake)

# Call the generation function with the provided arguments
generate_connection_attributes_registration("${SOURCE_FILE}" "${OUTPUT_FILE}")