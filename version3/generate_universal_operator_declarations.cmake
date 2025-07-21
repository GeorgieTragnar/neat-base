# Universal Operator Function Scanner
# Detects ALL functions in Operator namespace regardless of signature
# Generates both forward declarations and friend declarations

function(generate_universal_operator_declarations BASE_DIR DISCOVERED_TYPES_FILE FORWARDS_OUTPUT_FILE FRIENDS_OUTPUT_FILE)
    message(STATUS "Generating universal operator declarations from: ${BASE_DIR}")
    message(STATUS "Using discovered types from: ${DISCOVERED_TYPES_FILE}")
    message(STATUS "Forward declarations output: ${FORWARDS_OUTPUT_FILE}")
    message(STATUS "Friend declarations output: ${FRIENDS_OUTPUT_FILE}")
    
    # Load discovered types
    include(${DISCOVERED_TYPES_FILE})
    
    # Find all operator header files recursively
    file(GLOB_RECURSE OPERATOR_HEADERS "${BASE_DIR}/*.hpp")
    
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
                            string(FIND "${FINAL_OPERATOR_CONTENT}" "${CLASS_MATCH}" CLASS_MATCH_POS)
                            string(LENGTH "${CLASS_MATCH}" CLASS_MATCH_LEN)
                            math(EXPR SEARCH_START "${CLASS_MATCH_POS} + ${CLASS_MATCH_LEN}")
                            string(SUBSTRING "${FINAL_OPERATOR_CONTENT}" ${SEARCH_START} -1 AFTER_CLASS_CONTENT)
                            
                            # Extract first segment only
                            string(REGEX MATCH "([^;]*;)" FIRST_SEGMENT "${AFTER_CLASS_CONTENT}")
                            
                            if(FIRST_SEGMENT)
                                # Check if first segment is a param class declaration
                                string(REGEX MATCH "[ \t\n]*class[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t\n]*" FIRST_IS_CLASS "${FIRST_SEGMENT}")
                                
                                if(FIRST_IS_CLASS)
                                    # Extract second segment from remaining content after first semicolon
                                    string(REGEX REPLACE "^[^;]*;" "" REMAINING_CONTENT "${AFTER_CLASS_CONTENT}")
                                    string(REGEX MATCH "([^;]*;)" SECOND_SEGMENT "${REMAINING_CONTENT}")
                                    set(FUNCTION_CONTENT "${SECOND_SEGMENT}")
                                else()
                                    # Use first segment for function detection
                                    set(FUNCTION_CONTENT "${FIRST_SEGMENT}")
                                endif()
                                
                                # Apply type recognition to the selected segment
                                # Try template function first
                                string(REGEX MATCH "(template[ \t]*<[^>]*>)[ \t\n]*([A-Za-z_][A-Za-z0-9_:<>&*, \t\n]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)" 
                                       TEMPLATE_FUNCTION_MATCH "${FUNCTION_CONTENT}")
                                
                                if(TEMPLATE_FUNCTION_MATCH)
                                    set(TEMPLATE_DIRECTIVE "${CMAKE_MATCH_1}")
                                    set(RETURN_TYPE "${CMAKE_MATCH_2}")
                                    set(FUNCTION_NAME "${CMAKE_MATCH_3}")
                                    
                                    # Extract parameters by finding function name in full content and matching parentheses
                                    string(FIND "${FINAL_OPERATOR_CONTENT}" "${FUNCTION_NAME}" FUNC_NAME_POS)
                                    if(FUNC_NAME_POS GREATER -1)
                                        # Find opening parenthesis after function name
                                        string(LENGTH "${FUNCTION_NAME}" FUNC_NAME_LEN)
                                        math(EXPR SEARCH_START "${FUNC_NAME_POS} + ${FUNC_NAME_LEN}")
                                        string(SUBSTRING "${FINAL_OPERATOR_CONTENT}" ${SEARCH_START} -1 AFTER_FUNC_NAME)
                                        string(FIND "${AFTER_FUNC_NAME}" "(" FIRST_PAREN_POS)
                                        
                                        if(FIRST_PAREN_POS GREATER -1)
                                            # Count parentheses depth to find matching closing parenthesis
                                            # Also track angle brackets for template parameters
                                            math(EXPR PARAM_START "${FIRST_PAREN_POS} + 1")
                                            set(PAREN_DEPTH 1)
                                            set(ANGLE_DEPTH 0)
                                            set(CURRENT_POS ${PARAM_START})
                                            string(LENGTH "${AFTER_FUNC_NAME}" CONTENT_LENGTH)
                                            
                                            while(PAREN_DEPTH GREATER 0 AND CURRENT_POS LESS CONTENT_LENGTH)
                                                string(SUBSTRING "${AFTER_FUNC_NAME}" ${CURRENT_POS} 1 CHAR)
                                                if(CHAR STREQUAL "<")
                                                    math(EXPR ANGLE_DEPTH "${ANGLE_DEPTH} + 1")
                                                elseif(CHAR STREQUAL ">")
                                                    if(ANGLE_DEPTH GREATER 0)
                                                        math(EXPR ANGLE_DEPTH "${ANGLE_DEPTH} - 1")
                                                    endif()
                                                elseif(CHAR STREQUAL "(" AND ANGLE_DEPTH EQUAL 0)
                                                    # Only count parentheses when not inside template parameters
                                                    math(EXPR PAREN_DEPTH "${PAREN_DEPTH} + 1")
                                                elseif(CHAR STREQUAL ")" AND ANGLE_DEPTH EQUAL 0)
                                                    # Only count parentheses when not inside template parameters
                                                    math(EXPR PAREN_DEPTH "${PAREN_DEPTH} - 1")
                                                endif()
                                                math(EXPR CURRENT_POS "${CURRENT_POS} + 1")
                                            endwhile()
                                            
                                            if(PAREN_DEPTH EQUAL 0)
                                                math(EXPR PARAM_END "${CURRENT_POS} - 1")
                                                math(EXPR PARAM_LENGTH "${PARAM_END} - ${PARAM_START}")
                                                string(SUBSTRING "${AFTER_FUNC_NAME}" ${PARAM_START} ${PARAM_LENGTH} FUNCTION_PARAMS)
                                            else()
                                                set(FUNCTION_PARAMS "")
                                            endif()
                                        else()
                                            set(FUNCTION_PARAMS "")
                                        endif()
                                    else()
                                        set(FUNCTION_PARAMS "")
                                    endif()
                                    
                                    set(FUNCTION_FOUND TRUE)
                                    message(STATUS "  Found template operator function: ${FUNCTION_NAME}")
                                else()
                                    # Try regular function
                                    string(REGEX MATCH "([A-Za-z_][A-Za-z0-9_:<>&*, \t]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)" 
                                           REGULAR_FUNCTION_MATCH "${FUNCTION_CONTENT}")
                                    
                                    if(REGULAR_FUNCTION_MATCH)
                                        set(TEMPLATE_DIRECTIVE "")
                                        set(RETURN_TYPE "${CMAKE_MATCH_1}")
                                        set(FUNCTION_NAME "${CMAKE_MATCH_2}")
                                        
                                        # Extract parameters by finding function name in full content and matching parentheses
                                        string(FIND "${FINAL_OPERATOR_CONTENT}" "${FUNCTION_NAME}" FUNC_NAME_POS)
                                        if(FUNC_NAME_POS GREATER -1)
                                            # Find opening parenthesis after function name
                                            string(LENGTH "${FUNCTION_NAME}" FUNC_NAME_LEN)
                                            math(EXPR SEARCH_START "${FUNC_NAME_POS} + ${FUNC_NAME_LEN}")
                                            string(SUBSTRING "${FINAL_OPERATOR_CONTENT}" ${SEARCH_START} -1 AFTER_FUNC_NAME)
                                            string(FIND "${AFTER_FUNC_NAME}" "(" FIRST_PAREN_POS)
                                            
                                            if(FIRST_PAREN_POS GREATER -1)
                                                # Count parentheses depth to find matching closing parenthesis
                                                # Also track angle brackets for template parameters
                                                math(EXPR PARAM_START "${FIRST_PAREN_POS} + 1")
                                                set(PAREN_DEPTH 1)
                                                set(ANGLE_DEPTH 0)
                                                set(CURRENT_POS ${PARAM_START})
                                                string(LENGTH "${AFTER_FUNC_NAME}" CONTENT_LENGTH)
                                                
                                                while(PAREN_DEPTH GREATER 0 AND CURRENT_POS LESS CONTENT_LENGTH)
                                                    string(SUBSTRING "${AFTER_FUNC_NAME}" ${CURRENT_POS} 1 CHAR)
                                                    if(CHAR STREQUAL "<")
                                                        math(EXPR ANGLE_DEPTH "${ANGLE_DEPTH} + 1")
                                                    elseif(CHAR STREQUAL ">")
                                                        if(ANGLE_DEPTH GREATER 0)
                                                            math(EXPR ANGLE_DEPTH "${ANGLE_DEPTH} - 1")
                                                        endif()
                                                    elseif(CHAR STREQUAL "(" AND ANGLE_DEPTH EQUAL 0)
                                                        # Only count parentheses when not inside template parameters
                                                        math(EXPR PAREN_DEPTH "${PAREN_DEPTH} + 1")
                                                    elseif(CHAR STREQUAL ")" AND ANGLE_DEPTH EQUAL 0)
                                                        # Only count parentheses when not inside template parameters
                                                        math(EXPR PAREN_DEPTH "${PAREN_DEPTH} - 1")
                                                    endif()
                                                    math(EXPR CURRENT_POS "${CURRENT_POS} + 1")
                                                endwhile()
                                                
                                                if(PAREN_DEPTH EQUAL 0)
                                                    math(EXPR PARAM_END "${CURRENT_POS} - 1")
                                                    math(EXPR PARAM_LENGTH "${PARAM_END} - ${PARAM_START}")
                                                    string(SUBSTRING "${AFTER_FUNC_NAME}" ${PARAM_START} ${PARAM_LENGTH} FUNCTION_PARAMS)
                                                else()
                                                    set(FUNCTION_PARAMS "")
                                                endif()
                                            else()
                                                set(FUNCTION_PARAMS "")
                                            endif()
                                        else()
                                            set(FUNCTION_PARAMS "")
                                        endif()
                                        
                                        set(FUNCTION_FOUND TRUE)
                                        message(STATUS "  Found operator function: ${FUNCTION_NAME}")
                                    endif()
                                endif()
                            endif()
                        endif()
                        
                        # Pattern 2: Functions without Params class (fallback)
                        if(NOT FUNCTION_FOUND)
                            # Try template function first
                            string(REGEX MATCH "(template[ \t]*<[^>]*>)[ \t\n]*([A-Za-z_][A-Za-z0-9_:<>&*, \t]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\(([^)]*)\\)[ \t\n]*;" 
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
                                string(REGEX MATCH "([A-Za-z_][A-Za-z0-9_:<>&*, \t]+)[ \t\n]+([A-Za-z_][A-Za-z0-9_]*)[ \t\n]*\\((.*)\\)[ \t\n]*;" 
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
                            
                            # Create function signature for deduplication
                            set(FUNCTION_SIGNATURE "${CLEAN_RETURN_TYPE} ${CLEAN_FUNCTION_NAME}(${CLEAN_PARAMS})")
                            
                            # Check if we've already seen this function signature
                            set(ALREADY_EXISTS FALSE)
                            foreach(EXISTING_FUNCTION ${DISCOVERED_FUNCTIONS})
                                string(REPLACE "|" ";" EXISTING_PARTS "${EXISTING_FUNCTION}")
                                list(LENGTH EXISTING_PARTS EXISTING_PARTS_COUNT)
                                
                                if(EXISTING_PARTS_COUNT EQUAL 4)
                                    list(GET EXISTING_PARTS 1 EXISTING_RETURN_TYPE)
                                    list(GET EXISTING_PARTS 2 EXISTING_FUNCTION_NAME)
                                    list(GET EXISTING_PARTS 3 EXISTING_PARAMS)
                                    set(EXISTING_SIGNATURE "${EXISTING_RETURN_TYPE} ${EXISTING_FUNCTION_NAME}(${EXISTING_PARAMS})")
                                    
                                    if(FUNCTION_SIGNATURE STREQUAL EXISTING_SIGNATURE)
                                        set(ALREADY_EXISTS TRUE)
                                        break()
                                    endif()
                                endif()
                            endforeach()
                            
                            # Store the function only if it's not a duplicate and has valid data
                            if(NOT ALREADY_EXISTS AND CLEAN_FUNCTION_NAME AND CLEAN_RETURN_TYPE)
                                list(APPEND DISCOVERED_FUNCTIONS "${CLEAN_TEMPLATE_DIRECTIVE}|${CLEAN_RETURN_TYPE}|${CLEAN_FUNCTION_NAME}|${CLEAN_PARAMS}")
                            endif()
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
        else()
            # Skip malformed entries silently (expected behavior with current regex parsing)
            continue()
        endif()
            
        if(TEMPLATE_DIRECTIVE)
            string(APPEND FORWARDS_CONTENT "    ${TEMPLATE_DIRECTIVE}\n")
            string(APPEND FORWARDS_CONTENT "    ${RETURN_TYPE} ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
        else()
            string(APPEND FORWARDS_CONTENT "    ${RETURN_TYPE} ${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
        endif()
    endforeach()
    
    string(APPEND FORWARDS_CONTENT "}\n")
    
    # Postprocess forward declarations to fix extra )) from std::function parameters
    string(REPLACE "))" ")" FORWARDS_CONTENT "${FORWARDS_CONTENT}")
    
    # Postprocess to add missing template directives for functions with PopulationContainer<FitnessResultType>
    # Only add template directive if the line doesn't already have one
    string(REGEX REPLACE "\n([ \t]+[^<\n]*PopulationContainer<FitnessResultType>[^\n]*;)" 
           "\n    template<typename FitnessResultType> \\1" FORWARDS_CONTENT "${FORWARDS_CONTENT}")
    
    # Clean up any duplicate template directives that might have been created
    string(REGEX REPLACE "template<typename FitnessResultType>[ \t\n]+template<typename FitnessResultType>" 
           "template<typename FitnessResultType>" FORWARDS_CONTENT "${FORWARDS_CONTENT}")
    
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
        else()
            # Skip malformed entries silently (expected behavior with current regex parsing)
            continue()
        endif()
            
        if(TEMPLATE_DIRECTIVE)
            # For template functions, check if they contain FitnessResultType
            string(FIND "${FUNCTION_PARAMS}" "FitnessResultType" CONTAINS_FITNESS_TYPE)
            if(CONTAINS_FITNESS_TYPE GREATER -1)
                # Generate proper template friend declaration using the same parameter name
                string(APPEND FRIENDS_CONTENT "    template<typename FitnessResultType> friend ${RETURN_TYPE} Operator::${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            else()
                string(APPEND FRIENDS_CONTENT "    ${TEMPLATE_DIRECTIVE} friend ${RETURN_TYPE} Operator::${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
            endif()
        else()
            string(APPEND FRIENDS_CONTENT "    friend ${RETURN_TYPE} Operator::${FUNCTION_NAME}(${FUNCTION_PARAMS});\n")
        endif()
    endforeach()
    
    # Postprocess friend declarations to fix extra )) from std::function parameters
    string(REPLACE "))" ")" FRIENDS_CONTENT "${FRIENDS_CONTENT}")
    
    # Postprocess to add missing template directives for functions with PopulationContainer<FitnessResultType>
    # Only add template directive if the line doesn't already have one
    string(REGEX REPLACE "\n([ \t]+friend [^<\n]*PopulationContainer<FitnessResultType>[^\n]*;)" 
           "\n    template<typename FitnessResultType> \\1" FRIENDS_CONTENT "${FRIENDS_CONTENT}")
    
    # Clean up any duplicate template directives that might have been created
    string(REGEX REPLACE "template<typename FitnessResultType>[ \t\n]+template<typename FitnessResultType>" 
           "template<typename FitnessResultType>" FRIENDS_CONTENT "${FRIENDS_CONTENT}")
    
    # Write output files
    file(WRITE ${FORWARDS_OUTPUT_FILE} "${FORWARDS_CONTENT}")
    file(WRITE ${FRIENDS_OUTPUT_FILE} "${FRIENDS_CONTENT}")
    
    # Report results
    list(LENGTH DISCOVERED_FUNCTIONS FUNCTION_COUNT)
    message(STATUS "Generated declarations for ${FUNCTION_COUNT} operator functions")
    message(STATUS "Forward declarations saved to: ${FORWARDS_OUTPUT_FILE}")
    message(STATUS "Friend declarations saved to: ${FRIENDS_OUTPUT_FILE}")
endfunction()