# Function to generate NodeGeneAttributes field registration
function(generate_node_attributes_registration SOURCE_FILE OUTPUT_FILE)
    # Read the source file
    file(READ ${SOURCE_FILE} FILE_CONTENT)
    
    # Extract NodeGeneAttributes struct definition
    string(REGEX MATCH "struct[ \t\n]+NodeGeneAttributes[ \t\n]*{([^}]*)}" 
           STRUCT_MATCH "${FILE_CONTENT}")
    
    if(NOT STRUCT_MATCH)
        message(FATAL_ERROR "Could not find NodeGeneAttributes struct in ${SOURCE_FILE}")
    endif()
    
    # Get the content inside the struct
    set(STRUCT_CONTENT "${CMAKE_MATCH_1}")
    
    # Split into lines
    string(REGEX REPLACE ";" "\\\\;" STRUCT_CONTENT "${STRUCT_CONTENT}")
    string(REGEX REPLACE "\n" ";" STRUCT_LINES "${STRUCT_CONTENT}")
    
    # Initialize output
    set(REGISTRATION_CODE "")
    
    # Process each line
    foreach(LINE ${STRUCT_LINES})
        # Strip whitespace
        string(STRIP "${LINE}" LINE)
        
        # Skip empty lines and comments
        if(NOT LINE OR LINE MATCHES "^//")
            continue()
        endif()
        
        # Match member declarations (type name;)
        # This regex handles:
        # - Simple types: int32_t priority;
        # - Templated types: std::vector<int> values;
        # - Enum types: ActivationType activationType;
        string(REGEX MATCH "^([A-Za-z_][A-Za-z0-9_:<>]*(<[^>]+>)?)[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*;" 
               MEMBER_MATCH "${LINE}")
        
        if(MEMBER_MATCH)
            set(MEMBER_TYPE "${CMAKE_MATCH_1}")
            set(MEMBER_NAME "${CMAKE_MATCH_3}")
            
            # Generate registration line
            string(APPEND REGISTRATION_CODE 
                   "    register_field(create_field_descriptor(\"${MEMBER_NAME}\", &NodeGeneAttributes::${MEMBER_NAME}));\n")
        endif()
    endforeach()
    
    # Write output file
    file(WRITE ${OUTPUT_FILE} "${REGISTRATION_CODE}")
    
    # Report what was generated
    message(STATUS "Generated field registration for NodeGeneAttributes in ${OUTPUT_FILE}")
endfunction()