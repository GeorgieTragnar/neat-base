# Function to generate forward declarations for operators
function(generate_operator_forward_declarations OPERATOR_DIR OUTPUT_FILE)
    message(STATUS "Generating operator forwards from: ${OPERATOR_DIR}")
    message(STATUS "Output file: ${OUTPUT_FILE}")
    # Find all operator header files
    file(GLOB OPERATOR_HEADERS "${OPERATOR_DIR}/*.hpp")
    
    # Initialize output content
    set(FORWARD_DECLARATIONS "")
    set(FRIEND_DECLARATIONS "")
    set(NAMESPACE_DECLARATIONS "")
    set(PROCESSED_CLASSES "")
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
            string(REGEX MATCHALL "(template[ \t]*<[^>]*>[ \t\n]*)?Genome[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^)]*\\)" 
                   FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            # Also find void functions that take Genome& (phenotype operators)
            string(REGEX MATCHALL "void[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^)]*Genome[ \t]*&[^)]*\\)" 
                   VOID_FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            # Also find bool functions that take const Genome& (assertion operators)
            string(REGEX MATCHALL "bool[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\([^)]*const[ \t]+Genome[ \t]*&[^)]*\\)" 
                   BOOL_FUNCTION_MATCHES "${NAMESPACE_CONTENT}")
            
            foreach(FUNCTION_MATCH ${FUNCTION_MATCHES})
                # Extract function signature (including template part)
                string(REGEX MATCH "(template[ \t]*<[^>]*>[ \t\n]*)?Genome[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^)]*)\\)" 
                       SIGNATURE_MATCH "${FUNCTION_MATCH}")
                
                if(SIGNATURE_MATCH)
                    set(TEMPLATE_PART "${CMAKE_MATCH_1}")
                    set(FUNCTION_NAME "${CMAKE_MATCH_2}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_3}")
                    
                    # Extract template parameter names if this is a template function
                    set(TEMPLATE_PARAM_NAMES "")
                    if(TEMPLATE_PART)
                        string(REGEX MATCHALL "typename[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                               TEMPLATE_PARAM_MATCHES "${TEMPLATE_PART}")
                        foreach(TEMPLATE_PARAM_MATCH ${TEMPLATE_PARAM_MATCHES})
                            string(REGEX MATCH "typename[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                                   TEMPLATE_PARAM_EXTRACT "${TEMPLATE_PARAM_MATCH}")
                            if(TEMPLATE_PARAM_EXTRACT)
                                list(APPEND TEMPLATE_PARAM_NAMES "${CMAKE_MATCH_1}")
                            endif()
                        endforeach()
                    endif()
                    
                    # Find parameter class names, including namespaced ones (both const and non-const references)
                    string(REGEX MATCHALL "(const[ \t]+)?([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*&" 
                           PARAM_MATCHES "${FUNCTION_PARAMS}")
                    
                    foreach(PARAM_MATCH ${PARAM_MATCHES})
                        string(REGEX MATCH "(const[ \t]+)?([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*&" 
                               PARAM_CLASS_MATCH "${PARAM_MATCH}")
                        if(PARAM_CLASS_MATCH AND NOT CMAKE_MATCH_3 STREQUAL "Genome")
                            set(PARAM_NAMESPACE "${CMAKE_MATCH_2}")
                            set(PARAM_CLASS "${CMAKE_MATCH_3}")
                            
                            # Skip if this is a template parameter
                            list(FIND TEMPLATE_PARAM_NAMES "${PARAM_CLASS}" TEMPLATE_PARAM_INDEX)
                            if(TEMPLATE_PARAM_INDEX EQUAL -1)
                                # Create unique key for this class declaration
                                if(PARAM_NAMESPACE)
                                    string(REGEX REPLACE "::" "" CLEAN_NAMESPACE "${PARAM_NAMESPACE}")
                                    set(CLASS_KEY "${CLEAN_NAMESPACE}::${PARAM_CLASS}")
                                else()
                                    set(CLASS_KEY "::${PARAM_CLASS}")
                                endif()
                                
                                # Only add if not already processed
                                string(FIND "${PROCESSED_CLASSES}" "${CLASS_KEY};" CLASS_FOUND)
                                if(CLASS_FOUND EQUAL -1)
                                    string(APPEND PROCESSED_CLASSES "${CLASS_KEY};")
                                    
                                    # Add forward declaration for parameter class
                                    if(PARAM_NAMESPACE)
                                        string(APPEND NAMESPACE_DECLARATIONS 
                                               "namespace ${CLEAN_NAMESPACE} { class ${PARAM_CLASS}; }\n")
                                    else()
                                        string(APPEND FORWARD_DECLARATIONS 
                                               "\tclass ${PARAM_CLASS};\n")
                                    endif()
                                endif()
                            endif()
                        endif()
                    endforeach()
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "${FUNCTION_NAME}(${FUNCTION_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add function forward declaration (include template part if present)
                        if(TEMPLATE_PART)
                            string(APPEND FORWARD_DECLARATIONS 
                                   "\t${TEMPLATE_PART}Genome ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
                        else()
                            string(APPEND FORWARD_DECLARATIONS 
                                   "\tGenome ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
                        endif()
                    endif()
                endif()
            endforeach()
            
            # Process void functions that take Genome& (phenotype operators)
            foreach(VOID_FUNCTION_MATCH ${VOID_FUNCTION_MATCHES})
                # Extract function signature for void functions
                string(REGEX MATCH "void[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^)]*)\\)" 
                       VOID_SIGNATURE_MATCH "${VOID_FUNCTION_MATCH}")
                
                if(VOID_SIGNATURE_MATCH)
                    set(FUNCTION_NAME "${CMAKE_MATCH_1}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_2}")
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "${FUNCTION_NAME}(${FUNCTION_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add void function forward declaration
                        string(APPEND FORWARD_DECLARATIONS 
                               "\tvoid ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
                    endif()
                endif()
            endforeach()
            
            # Process bool functions that take const Genome& (assertion operators)
            foreach(BOOL_FUNCTION_MATCH ${BOOL_FUNCTION_MATCHES})
                # Extract function signature for bool functions
                string(REGEX MATCH "bool[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(([^)]*)\\)" 
                       BOOL_SIGNATURE_MATCH "${BOOL_FUNCTION_MATCH}")
                
                if(BOOL_SIGNATURE_MATCH)
                    set(FUNCTION_NAME "${CMAKE_MATCH_1}")
                    set(FUNCTION_PARAMS "${CMAKE_MATCH_2}")
                    
                    # Find parameter class names, including namespaced ones (both const and non-const references)
                    string(REGEX MATCHALL "(const[ \t]+)?([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*&" 
                           PARAM_MATCHES "${FUNCTION_PARAMS}")
                    
                    foreach(PARAM_MATCH ${PARAM_MATCHES})
                        string(REGEX MATCH "(const[ \t]+)?([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*&" 
                               PARAM_CLASS_MATCH "${PARAM_MATCH}")
                        if(PARAM_CLASS_MATCH AND NOT CMAKE_MATCH_3 STREQUAL "Genome")
                            set(PARAM_NAMESPACE "${CMAKE_MATCH_2}")
                            set(PARAM_CLASS "${CMAKE_MATCH_3}")
                            
                            # Create unique key for this class declaration
                            if(PARAM_NAMESPACE)
                                string(REGEX REPLACE "::" "" CLEAN_NAMESPACE "${PARAM_NAMESPACE}")
                                set(CLASS_KEY "${CLEAN_NAMESPACE}::${PARAM_CLASS}")
                            else()
                                set(CLASS_KEY "::${PARAM_CLASS}")
                            endif()
                            
                            # Only add if not already processed
                            string(FIND "${PROCESSED_CLASSES}" "${CLASS_KEY};" CLASS_FOUND)
                            if(CLASS_FOUND EQUAL -1)
                                string(APPEND PROCESSED_CLASSES "${CLASS_KEY};")
                                
                                # Add forward declaration for parameter class
                                if(PARAM_NAMESPACE)
                                    string(APPEND NAMESPACE_DECLARATIONS 
                                           "namespace ${CLEAN_NAMESPACE} { class ${PARAM_CLASS}; }\n")
                                else()
                                    string(APPEND FORWARD_DECLARATIONS 
                                           "\tclass ${PARAM_CLASS};\n")
                                endif()
                            endif()
                        endif()
                    endforeach()
                    
                    # Create unique key for this function
                    set(FUNCTION_KEY "${FUNCTION_NAME}(${FUNCTION_PARAMS})")
                    
                    # Only add if not already processed
                    string(FIND "${PROCESSED_FUNCTIONS}" "${FUNCTION_KEY};" FUNCTION_FOUND)
                    if(FUNCTION_FOUND EQUAL -1)
                        string(APPEND PROCESSED_FUNCTIONS "${FUNCTION_KEY};")
                        # Add bool function forward declaration
                        string(APPEND FORWARD_DECLARATIONS 
                               "\tbool ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
                    endif()
                endif()
            endforeach()
        endif()
    endforeach()
    
    # Generate complete include content
    set(INCLUDE_CONTENT "// Auto-generated forward declarations for operators\n")
    string(APPEND INCLUDE_CONTENT "// This file is generated by CMake - do not edit manually\n\n")
    string(APPEND INCLUDE_CONTENT "class Genome;\n")
    string(APPEND INCLUDE_CONTENT "class HistoryTracker;\n\n")
    # Add namespace declarations outside of Operator namespace
    if(NAMESPACE_DECLARATIONS)
        string(APPEND INCLUDE_CONTENT "${NAMESPACE_DECLARATIONS}\n")
    endif()
    string(APPEND INCLUDE_CONTENT "namespace Operator {\n")
    string(APPEND INCLUDE_CONTENT "${FORWARD_DECLARATIONS}")
    string(APPEND INCLUDE_CONTENT "}\n")
    
    # Write output file
    file(WRITE ${OUTPUT_FILE} "${INCLUDE_CONTENT}")
    
    # Report what was generated
    message(STATUS "Generated operator forward declarations in ${OUTPUT_FILE}")
endfunction()