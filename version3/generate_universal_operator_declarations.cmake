# Universal Operator Function Scanner
# Detects ALL functions in Operator namespace regardless of signature
# Generates both forward declarations and friend declarations

function(generate_universal_operator_declarations OPERATOR_DIR DISCOVERED_TYPES_FILE FORWARDS_OUTPUT_FILE FRIENDS_OUTPUT_FILE)
    message(STATUS "Generating universal operator declarations from: ${OPERATOR_DIR}")
    message(STATUS "Using discovered types from: ${DISCOVERED_TYPES_FILE}")
    message(STATUS "Forward declarations output: ${FORWARDS_OUTPUT_FILE}")
    message(STATUS "Friend declarations output: ${FRIENDS_OUTPUT_FILE}")
    
    # Load discovered types
    include(${DISCOVERED_TYPES_FILE})
    
    # Find all operator header files
    file(GLOB OPERATOR_HEADERS "${OPERATOR_DIR}/*.hpp")
    
    # Initialize collections
    set(DISCOVERED_FUNCTIONS "")
    
    # Process each operator header
    foreach(HEADER_FILE ${OPERATOR_HEADERS})
        get_filename_component(HEADER_NAME ${HEADER_FILE} NAME)
        message(STATUS "Scanning operator file: ${HEADER_NAME}")
        
        file(READ ${HEADER_FILE} FILE_CONTENT)
        
        # Check if file contains namespace Operator
        string(FIND "${FILE_CONTENT}" "namespace Operator" NAMESPACE_POS)
        
        if(NAMESPACE_POS GREATER -1)
            # Remove anonymous namespace blocks to avoid false positives
            set(CLEAN_CONTENT "${FILE_CONTENT}")
            
            # Find and remove all anonymous namespace blocks
            string(REGEX REPLACE "namespace[ \t]*\\{[^}]*\\}" "" CLEAN_CONTENT "${CLEAN_CONTENT}")
            
            # Extract only the Operator namespace content
            string(FIND "${CLEAN_CONTENT}" "namespace Operator" OP_NAMESPACE_START)
            if(OP_NAMESPACE_START GREATER -1)
                string(SUBSTRING "${CLEAN_CONTENT}" ${OP_NAMESPACE_START} -1 OP_NAMESPACE_SECTION)
                string(FIND "${OP_NAMESPACE_SECTION}" "{" BRACE_START)
                if(BRACE_START GREATER -1)
                    math(EXPR CONTENT_START "${BRACE_START} + 1")
                    string(SUBSTRING "${OP_NAMESPACE_SECTION}" ${CONTENT_START} -1 OPERATOR_CONTENT)
                    
                    # Find the matching closing brace
                    set(BRACE_COUNT 1)
                    set(POS 0)
                    string(LENGTH "${OPERATOR_CONTENT}" CONTENT_LENGTH)
                    
                    while(BRACE_COUNT GREATER 0 AND POS LESS CONTENT_LENGTH)
                        string(SUBSTRING "${OPERATOR_CONTENT}" ${POS} 1 CHAR)
                        if(CHAR STREQUAL "{")
                            math(EXPR BRACE_COUNT "${BRACE_COUNT} + 1")
                        elseif(CHAR STREQUAL "}")
                            math(EXPR BRACE_COUNT "${BRACE_COUNT} - 1")
                        endif()
                        math(EXPR POS "${POS} + 1")
                    endwhile()
                    
                    if(BRACE_COUNT EQUAL 0)
                        math(EXPR END_POS "${POS} - 1")
                        string(SUBSTRING "${OPERATOR_CONTENT}" 0 ${END_POS} FINAL_OPERATOR_CONTENT)
                        
                        # Enhanced pattern-based detection: Handle multiple patterns
                        set(FUNCTION_FOUND FALSE)
                        
                        # Pattern 1: Functions with Params class forward declaration
                        string(REGEX MATCH "class[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]*;" CLASS_MATCH "${FINAL_OPERATOR_CONTENT}")
                        
                        if(CLASS_MATCH)
                            string(FIND "${FINAL_OPERATOR_CONTENT}" "${CLASS_MATCH}" CLASS_END_POS)
                            math(EXPR SEARCH_START "${CLASS_END_POS} + 1")
                            string(SUBSTRING "${FINAL_OPERATOR_CONTENT}" ${SEARCH_START} -1 AFTER_CLASS_CONTENT)
                            
                            # Try template function first
                            string(REGEX MATCH "(template[ \t]*<[^>]*>)[ \t\n]*([A-Za-z_][A-Za-z0-9_:<>&*]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\(([^)]*)\\)[ \t\n]*;" 
                                   TEMPLATE_FUNCTION_MATCH "${AFTER_CLASS_CONTENT}")
                            
                            if(TEMPLATE_FUNCTION_MATCH)
                                set(TEMPLATE_DIRECTIVE "${CMAKE_MATCH_1}")
                                set(RETURN_TYPE "${CMAKE_MATCH_2}")
                                set(FUNCTION_NAME "${CMAKE_MATCH_3}")
                                set(FUNCTION_PARAMS "${CMAKE_MATCH_4}")
                                set(FUNCTION_FOUND TRUE)
                                message(STATUS "  Found template operator function: ${FUNCTION_NAME}")
                            else()
                                # Try regular function
                                string(REGEX MATCH "([A-Za-z_][A-Za-z0-9_:<>&*]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\(([^)]*)\\)[ \t\n]*;" 
                                       REGULAR_FUNCTION_MATCH "${AFTER_CLASS_CONTENT}")
                                
                                if(REGULAR_FUNCTION_MATCH)
                                    set(TEMPLATE_DIRECTIVE "")
                                    set(RETURN_TYPE "${CMAKE_MATCH_1}")
                                    set(FUNCTION_NAME "${CMAKE_MATCH_2}")
                                    set(FUNCTION_PARAMS "${CMAKE_MATCH_3}")
                                    set(FUNCTION_FOUND TRUE)
                                    message(STATUS "  Found operator function: ${FUNCTION_NAME}")
                                endif()
                            endif()
                        endif()
                        
                        # Pattern 2: Functions without Params class (fallback)
                        if(NOT FUNCTION_FOUND)
                            # Try template function first
                            string(REGEX MATCH "(template[ \t]*<[^>]*>)[ \t\n]*([A-Za-z_][A-Za-z0-9_:<>&*]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\(([^)]*)\\)[ \t\n]*;" 
                                   TEMPLATE_FUNCTION_MATCH "${FINAL_OPERATOR_CONTENT}")
                            
                            if(TEMPLATE_FUNCTION_MATCH)
                                set(TEMPLATE_DIRECTIVE "${CMAKE_MATCH_1}")
                                set(RETURN_TYPE "${CMAKE_MATCH_2}")
                                set(FUNCTION_NAME "${CMAKE_MATCH_3}")
                                set(FUNCTION_PARAMS "${CMAKE_MATCH_4}")
                                set(FUNCTION_FOUND TRUE)
                                message(STATUS "  Found template operator function (no params class): ${FUNCTION_NAME}")
                            else()
                                # Try regular function
                                string(REGEX MATCH "([A-Za-z_][A-Za-z0-9_:<>&*]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\(([^)]*)\\)[ \t\n]*;" 
                                       REGULAR_FUNCTION_MATCH "${FINAL_OPERATOR_CONTENT}")
                                
                                if(REGULAR_FUNCTION_MATCH)
                                    set(TEMPLATE_DIRECTIVE "")
                                    set(RETURN_TYPE "${CMAKE_MATCH_1}")
                                    set(FUNCTION_NAME "${CMAKE_MATCH_2}")
                                    set(FUNCTION_PARAMS "${CMAKE_MATCH_3}")
                                    set(FUNCTION_FOUND TRUE)
                                    message(STATUS "  Found operator function (no params class): ${FUNCTION_NAME}")
                                endif()
                            endif()
                        endif()
                        
                        if(FUNCTION_FOUND)
                            # Clean up whitespace and normalize parameters
                            string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_TEMPLATE_DIRECTIVE "${TEMPLATE_DIRECTIVE}")
                            string(STRIP "${CLEAN_TEMPLATE_DIRECTIVE}" CLEAN_TEMPLATE_DIRECTIVE)
                            string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_RETURN_TYPE "${RETURN_TYPE}")
                            string(STRIP "${CLEAN_RETURN_TYPE}" CLEAN_RETURN_TYPE)
                            string(STRIP "${FUNCTION_NAME}" CLEAN_FUNCTION_NAME)
                            string(REGEX REPLACE "[ \t\n\r]+" " " CLEAN_PARAMS "${FUNCTION_PARAMS}")
                            string(STRIP "${CLEAN_PARAMS}" CLEAN_PARAMS)
                            
                            # Store the function found with template directive
                            list(APPEND DISCOVERED_FUNCTIONS "${CLEAN_TEMPLATE_DIRECTIVE}|${CLEAN_RETURN_TYPE}|${CLEAN_FUNCTION_NAME}|${CLEAN_PARAMS}")
                        else()
                            message(STATUS "  No function found matching any expected pattern")
                        endif()
                    else()
                        message(STATUS "  Could not parse namespace content properly")
                    endif()
                else()
                    message(STATUS "  Could not find opening brace for namespace")
                endif()
            else()
                message(STATUS "  Could not find namespace Operator after cleaning")
            endif()
        endif()
    endforeach()
    
    # Generate forward declarations file
    set(FORWARDS_CONTENT "// Auto-generated operator forward declarations\n")
    string(APPEND FORWARDS_CONTENT "// This file is generated by CMake - do not edit manually\n\n")
    string(APPEND FORWARDS_CONTENT "namespace Operator {\n")
    
    foreach(FUNCTION_INFO ${DISCOVERED_FUNCTIONS})
        string(REPLACE "|" ";" FUNC_PARTS "${FUNCTION_INFO}")
        list(LENGTH FUNC_PARTS PARTS_COUNT)
        
        if(PARTS_COUNT EQUAL 4)
            list(GET FUNC_PARTS 0 TEMPLATE_DIRECTIVE)
            list(GET FUNC_PARTS 1 RETURN_TYPE)
            list(GET FUNC_PARTS 2 FUNCTION_NAME)
            list(GET FUNC_PARTS 3 FUNCTION_PARAMS)
            
            if(TEMPLATE_DIRECTIVE)
                string(APPEND FORWARDS_CONTENT "    ${TEMPLATE_DIRECTIVE}\n")
                string(APPEND FORWARDS_CONTENT "    ${RETURN_TYPE} ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            else()
                string(APPEND FORWARDS_CONTENT "    ${RETURN_TYPE} ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            endif()
        endif()
    endforeach()
    
    string(APPEND FORWARDS_CONTENT "}\n")
    
    # Generate friend declarations file
    set(FRIENDS_CONTENT "// Auto-generated operator friend declarations\n")
    string(APPEND FRIENDS_CONTENT "// This file is generated by CMake - do not edit manually\n")
    string(APPEND FRIENDS_CONTENT "// Include this in the protected section of data structure classes\n\n")
    
    foreach(FUNCTION_INFO ${DISCOVERED_FUNCTIONS})
        string(REPLACE "|" ";" FUNC_PARTS "${FUNCTION_INFO}")
        list(LENGTH FUNC_PARTS PARTS_COUNT)
        
        if(PARTS_COUNT EQUAL 4)
            list(GET FUNC_PARTS 0 TEMPLATE_DIRECTIVE)
            list(GET FUNC_PARTS 1 RETURN_TYPE)
            list(GET FUNC_PARTS 2 FUNCTION_NAME)
            list(GET FUNC_PARTS 3 FUNCTION_PARAMS)
            
            if(TEMPLATE_DIRECTIVE)
                string(APPEND FRIENDS_CONTENT "    ${TEMPLATE_DIRECTIVE} friend ${RETURN_TYPE} Operator::${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            else()
                string(APPEND FRIENDS_CONTENT "    friend ${RETURN_TYPE} Operator::${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            endif()
        endif()
    endforeach()
    
    # Write output files
    file(WRITE ${FORWARDS_OUTPUT_FILE} "${FORWARDS_CONTENT}")
    file(WRITE ${FRIENDS_OUTPUT_FILE} "${FRIENDS_CONTENT}")
    
    # Report results
    list(LENGTH DISCOVERED_FUNCTIONS FUNCTION_COUNT)
    message(STATUS "Generated declarations for ${FUNCTION_COUNT} operator functions")
    message(STATUS "Forward declarations saved to: ${FORWARDS_OUTPUT_FILE}")
    message(STATUS "Friend declarations saved to: ${FRIENDS_OUTPUT_FILE}")
endfunction()