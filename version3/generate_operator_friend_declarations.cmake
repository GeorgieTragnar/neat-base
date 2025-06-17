# Function to generate friend declarations for operators
function(generate_operator_friend_declarations OPERATOR_DIR OUTPUT_FILE)
    message(STATUS "Generating operator friend declarations from: ${OPERATOR_DIR}")
    message(STATUS "Output file: ${OUTPUT_FILE}")
    
    # Find all operator header files
    file(GLOB OPERATOR_HEADERS "${OPERATOR_DIR}/*.hpp")
    
    # Initialize output content
    set(FRIEND_DECLARATIONS "")
    set(PROCESSED_FUNCTIONS "")
    
    # Process each operator header
    foreach(HEADER_FILE ${OPERATOR_HEADERS})
        # Read the header file
        file(READ ${HEADER_FILE} FILE_CONTENT)
        
        # Extract namespace Operator content
        string(REGEX MATCH "namespace[ \t]+Operator[ \t]*{([^}]*)}" 
               NAMESPACE_MATCH "${FILE_CONTENT}")
        
        if(NAMESPACE_MATCH)
            set(NAMESPACE_CONTENT "${CMAKE_MATCH_1}")
            
            # Find function declarations that return Genome (including template functions)
            # Updated regex to handle multiline function declarations and templates
            string(REGEX MATCHALL "(template[ \t]*<[^>]*>[ \t\n]*)?Genome[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^{;]*\\)" 
                   FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            # Also find void functions that take Genome& (phenotype operators)
            string(REGEX MATCHALL "void[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^{;]*Genome[ \t]*&[^{;]*\\)" 
                   VOID_FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            # Also find bool functions that take const Genome& (assertion operators)
            string(REGEX MATCHALL "bool[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^{;]*const[ \t]+Genome[ \t]*&[^{;]*\\)" 
                   BOOL_FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            foreach(FUNCTION_MATCH ${FUNCTION_MATCHES})
                # Extract function signature with better multiline handling (including template part)
                string(REGEX MATCH "(template[ \t]*<[^>]*>[ \t\n]*)?Genome[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^{;]*)\\)" 
                       SIGNATURE_MATCH "${FUNCTION_MATCH}")
                
                if(SIGNATURE_MATCH)
                    set(TEMPLATE_PART "${CMAKE_MATCH_1}")
                    set(FUNCTION_NAME "${CMAKE_MATCH_2}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_3}")
                    
                    # Clean up parameters (remove extra whitespace and newlines)
                    string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_PARAMS "${FUNCTION_PARAMS}")
                    string(STRIP "${CLEAN_PARAMS}" CLEAN_PARAMS)
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "Genome ${FUNCTION_NAME}(${CLEAN_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add friend declaration (fully qualified with namespace, include template part if present)
                        if(TEMPLATE_PART)
                            # Clean up template part
                            string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_TEMPLATE "${TEMPLATE_PART}")
                            string(STRIP "${CLEAN_TEMPLATE}" CLEAN_TEMPLATE)
                            string(APPEND FRIEND_DECLARATIONS 
                                   "\t${CLEAN_TEMPLATE} friend Genome Operator::${FUNCTION_NAME}(${CLEAN_PARAMS});\n")
                        else()
                            string(APPEND FRIEND_DECLARATIONS 
                                   "\tfriend Genome Operator::${FUNCTION_NAME}(${CLEAN_PARAMS});\n")
                        endif()
                    endif()
                endif()
            endforeach()
            
            # Process void functions that take Genome& (phenotype operators)
            foreach(VOID_FUNCTION_MATCH ${VOID_FUNCTION_MATCHES})
                # Extract function signature for void functions
                string(REGEX MATCH "void[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^{;]*)\\)" 
                       VOID_SIGNATURE_MATCH "${VOID_FUNCTION_MATCH}")
                
                if(VOID_SIGNATURE_MATCH)
                    set(FUNCTION_NAME "${CMAKE_MATCH_1}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_2}")
                    
                    # Clean up parameters (remove extra whitespace and newlines)
                    string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_PARAMS "${FUNCTION_PARAMS}")
                    string(STRIP "${CLEAN_PARAMS}" CLEAN_PARAMS)
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "void ${FUNCTION_NAME}(${CLEAN_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add friend declaration for void function
                        string(APPEND FRIEND_DECLARATIONS 
                               "\tfriend void Operator::${FUNCTION_NAME}(${CLEAN_PARAMS});\n")
                    endif()
                endif()
            endforeach()
            
            # Process bool functions that take const Genome& (assertion operators)
            foreach(BOOL_FUNCTION_MATCH ${BOOL_FUNCTION_MATCHES})
                # Extract function signature for bool functions
                string(REGEX MATCH "bool[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^{;]*)\\)" 
                       BOOL_SIGNATURE_MATCH "${BOOL_FUNCTION_MATCH}")
                
                if(BOOL_SIGNATURE_MATCH)
                    set(FUNCTION_NAME "${CMAKE_MATCH_1}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_2}")
                    
                    # Clean up parameters (remove extra whitespace and newlines)
                    string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_PARAMS "${FUNCTION_PARAMS}")
                    string(STRIP "${CLEAN_PARAMS}" CLEAN_PARAMS)
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "bool ${FUNCTION_NAME}(${CLEAN_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add friend declaration for bool function
                        string(APPEND FRIEND_DECLARATIONS 
                               "\tfriend bool Operator::${FUNCTION_NAME}(${CLEAN_PARAMS});\n")
                    endif()
                endif()
            endforeach()
        endif()
    endforeach()
    
    # Generate complete include content
    set(INCLUDE_CONTENT "// Auto-generated friend declarations for operators\n")
    string(APPEND INCLUDE_CONTENT "// This file is generated by CMake - do not edit manually\n")
    string(APPEND INCLUDE_CONTENT "// Include this in the protected section of data structure classes\n\n")
    string(APPEND INCLUDE_CONTENT "${FRIEND_DECLARATIONS}")
    
    # Write output file
    file(WRITE ${OUTPUT_FILE} "${INCLUDE_CONTENT}")
    
    # Report what was generated
    message(STATUS "Generated operator friend declarations in ${OUTPUT_FILE}")
endfunction()