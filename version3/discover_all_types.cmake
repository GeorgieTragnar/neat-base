# Universal Type Discovery Scanner
# Scans all library code and discovers all class/struct types dynamically
# Excludes STL types and handles namespaces properly

function(discover_all_types SOURCE_DIR OUTPUT_FILE)
    message(STATUS "Discovering all types in: ${SOURCE_DIR}")
    message(STATUS "Output file: ${OUTPUT_FILE}")
    
    # Find all header files in all directories under version3
    file(GLOB_RECURSE ALL_HEADERS 
         "${SOURCE_DIR}/*.hpp")
    
    # Initialize discovery containers
    set(ALL_DISCOVERED_TYPES "")
    set(GLOBAL_TYPES "")
    set(POPULATION_TYPES "")
    set(OPERATOR_TYPES "")
    set(TEMPLATE_TYPES "")
    
    # Process each header file
    foreach(HEADER_FILE ${ALL_HEADERS})
        message(STATUS "Scanning: ${HEADER_FILE}")
        file(READ ${HEADER_FILE} FILE_CONTENT)
        
        # Discover types from multiple sources:
        
        # 1. Class/struct declarations
        string(REGEX MATCHALL "(class|struct)[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
               CLASS_DECLARATIONS "${FILE_CONTENT}")
        
        # 2. Template class declarations  
        string(REGEX MATCHALL "template[ \t]*<[^>]*>[ \t\n]*(class|struct)[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
               TEMPLATE_CLASS_DECLARATIONS "${FILE_CONTENT}")
        
        # 3. Function parameter types (reference and pointer)
        string(REGEX MATCHALL "([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*[&*]" 
               PARAM_TYPE_REFERENCES "${FILE_CONTENT}")
        
        # 4. Function return types
        string(REGEX MATCHALL "^[ \t]*([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\\(" 
               RETURN_TYPES "${FILE_CONTENT}")
        
        # 5. Template instantiations (e.g., PopulationContainer<FitnessResultType>)
        string(REGEX MATCHALL "([A-Za-z_][A-Za-z0-9_]*)<([A-Za-z_][A-Za-z0-9_]*)>" 
               TEMPLATE_INSTANTIATIONS "${FILE_CONTENT}")
        
        # Process class declarations
        foreach(CLASS_DECL ${CLASS_DECLARATIONS})
            string(REGEX MATCH "(class|struct)[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                   CLASS_MATCH "${CLASS_DECL}")
            if(CLASS_MATCH)
                set(CLASS_NAME "${CMAKE_MATCH_2}")
                
                # Skip STL types and invalid identifiers
                if(NOT CLASS_NAME MATCHES "^(string|vector|map|unordered_map|set|unordered_set|list|queue|stack|shared_ptr|unique_ptr|auto|enabled|template|typename|namespace|const|static|inline|virtual|override|public|private|protected)$")
                    list(APPEND ALL_DISCOVERED_TYPES "${CLASS_NAME}")
                endif()
            endif()
        endforeach()
        
        # Process template class declarations
        foreach(TEMPLATE_DECL ${TEMPLATE_CLASS_DECLARATIONS})
            string(REGEX MATCH "template[ \t]*<([^>]*)>[ \t\n]*(class|struct)[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                   TEMPLATE_MATCH "${TEMPLATE_DECL}")
            if(TEMPLATE_MATCH)
                set(TEMPLATE_PARAMS "${CMAKE_MATCH_1}")
                set(CLASS_NAME "${CMAKE_MATCH_3}")
                
                # Extract template parameter names
                string(REGEX MATCHALL "typename[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                       TEMPLATE_PARAM_MATCHES "${TEMPLATE_PARAMS}")
                
                set(TEMPLATE_PARAM_LIST "")
                foreach(TEMPLATE_PARAM_MATCH ${TEMPLATE_PARAM_MATCHES})
                    string(REGEX MATCH "typename[ \t]+([A-Za-z_][A-Za-z0-9_]*)" 
                           PARAM_EXTRACT "${TEMPLATE_PARAM_MATCH}")
                    if(PARAM_EXTRACT)
                        if(TEMPLATE_PARAM_LIST)
                            set(TEMPLATE_PARAM_LIST "${TEMPLATE_PARAM_LIST}, ${CMAKE_MATCH_1}")
                        else()
                            set(TEMPLATE_PARAM_LIST "${CMAKE_MATCH_1}")
                        endif()
                    endif()
                endforeach()
                
                # Create template type entry
                if(TEMPLATE_PARAM_LIST)
                    list(APPEND TEMPLATE_TYPES "template<typename ${TEMPLATE_PARAM_LIST}> class ${CLASS_NAME}")
                endif()
            endif()
        endforeach()
        
        # Process parameter type references
        foreach(PARAM_REF ${PARAM_TYPE_REFERENCES})
            string(REGEX MATCH "([A-Za-z_][A-Za-z0-9_]*::)?([A-Za-z_][A-Za-z0-9_]*)[ \t]*[&*]" 
                   PARAM_MATCH "${PARAM_REF}")
            if(PARAM_MATCH)
                set(PARAM_NAMESPACE "${CMAKE_MATCH_1}")
                set(PARAM_CLASS "${CMAKE_MATCH_2}")
                
                # Skip STL, primitive types, and keywords
                if(NOT PARAM_CLASS MATCHES "^(std|string|vector|map|unordered_map|set|unordered_set|list|queue|stack|shared_ptr|unique_ptr|int|uint32_t|size_t|bool|double|float|char|auto|enabled|template|typename|namespace|const|static|inline|virtual|override|public|private|protected|speciesSize|void)$" AND
                   NOT PARAM_NAMESPACE STREQUAL "std::")
                    list(APPEND ALL_DISCOVERED_TYPES "${PARAM_CLASS}")
                endif()
            endif()
        endforeach()
    endforeach()
    
    # Remove duplicates
    list(REMOVE_DUPLICATES ALL_DISCOVERED_TYPES)
    if(TEMPLATE_TYPES)
        list(REMOVE_DUPLICATES TEMPLATE_TYPES)
    endif()
    
    # Get explicit file lists for categorization
    file(GLOB_RECURSE ALL_VERSION3_HEADERS "${SOURCE_DIR}/*.hpp")
    
    # Categorize types by proper C++ scope parsing
    foreach(TYPE_NAME ${ALL_DISCOVERED_TYPES})
        set(TYPE_FOUND FALSE)
        
        # Check if type is defined in any version3 file
        foreach(HEADER_FILE ${ALL_VERSION3_HEADERS})
            file(READ "${HEADER_FILE}" FILE_CONTENT)
            
            # Look for the actual definition (not just any occurrence)
            string(REGEX MATCH "(^|[\n])[ \t]*(class|struct|enum[ \t]+class)[ \t]+${TYPE_NAME}[ \t]*[{;]" DEF_MATCH "${FILE_CONTENT}")
            string(REGEX MATCH "namespace[ \t]+Population[ \t]*\\{[^}]*${DEF_MATCH}" POP_NAMESPACE_MATCH "${FILE_CONTENT}")
            string(REGEX MATCH "namespace[ \t]+Operator[ \t]*\\{[^}]*${DEF_MATCH}" OP_NAMESPACE_MATCH "${FILE_CONTENT}")
            
            if(DEF_MATCH)
                # Extract the type keyword (class, struct, enum class)
                string(REGEX MATCH "(class|struct|enum[ \t]+class)" TYPE_KEYWORD "${DEF_MATCH}")
                
                # Check if this is a nested class by looking for class definition context
                # Find the position of this type definition
                string(FIND "${FILE_CONTENT}" "${DEF_MATCH}" DEF_POS)
                
                # Look backwards from the definition to see if we're inside another class
                string(SUBSTRING "${FILE_CONTENT}" 0 ${DEF_POS} CONTENT_BEFORE)
                
                # Count unmatched opening braces from class definitions
                string(REGEX MATCHALL "(^|[\n])[ \t]*(class|struct)[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]*\\{" CLASS_OPENS "${CONTENT_BEFORE}")
                string(REGEX MATCHALL "\\}" BRACE_CLOSES "${CONTENT_BEFORE}")
                
                list(LENGTH CLASS_OPENS OPEN_COUNT)
                list(LENGTH BRACE_CLOSES CLOSE_COUNT)
                
                # If more class opens than total closes, we're nested
                if(OPEN_COUNT GREATER CLOSE_COUNT)
                    message(STATUS "  Skipping nested class ${TYPE_NAME} (cannot be forward declared)")
                    # Don't add nested types to any category
                else()
                    if(POP_NAMESPACE_MATCH)
                        list(APPEND POPULATION_TYPES "${TYPE_NAME}")
                        message(STATUS "  Found ${TYPE_KEYWORD} ${TYPE_NAME} in Population namespace")
                    elseif(OP_NAMESPACE_MATCH)
                        list(APPEND OPERATOR_TYPES "${TYPE_NAME}")
                        message(STATUS "  Found ${TYPE_KEYWORD} ${TYPE_NAME} in Operator namespace")
                    else()
                        # Check if it's an enum class (needs special handling)
                        if(TYPE_KEYWORD MATCHES "enum[ \t]+class")
                            list(APPEND GLOBAL_TYPES "${TYPE_NAME}:ENUM")
                            message(STATUS "  Found enum class ${TYPE_NAME} in global scope")
                        else()
                            list(APPEND GLOBAL_TYPES "${TYPE_NAME}")
                            message(STATUS "  Found ${TYPE_KEYWORD} ${TYPE_NAME} in global scope")
                        endif()
                    endif()
                endif()
                set(TYPE_FOUND TRUE)
                break()
            endif()
        endforeach()
        
        if(NOT TYPE_FOUND)
            message(STATUS "  Warning: Could not determine namespace for type ${TYPE_NAME}")
        endif()
    endforeach()
    
    # Generate output file with discovered types
    set(OUTPUT_CONTENT "# Auto-generated type discovery results\n")
    string(APPEND OUTPUT_CONTENT "# This file contains all discovered types categorized by namespace\n\n")
    
    # Global types
    if(GLOBAL_TYPES)
        string(APPEND OUTPUT_CONTENT "set(DISCOVERED_GLOBAL_TYPES \"")
        string(REPLACE ";" ";" GLOBAL_TYPES_STR "${GLOBAL_TYPES}")
        string(APPEND OUTPUT_CONTENT "${GLOBAL_TYPES_STR}\")\n\n")
    endif()
    
    # Population namespace types  
    if(POPULATION_TYPES)
        string(APPEND OUTPUT_CONTENT "set(DISCOVERED_POPULATION_TYPES \"")
        string(REPLACE ";" ";" POPULATION_TYPES_STR "${POPULATION_TYPES}")
        string(APPEND OUTPUT_CONTENT "${POPULATION_TYPES_STR}\")\n\n")
    endif()
    
    # Operator namespace types
    if(OPERATOR_TYPES)
        string(APPEND OUTPUT_CONTENT "set(DISCOVERED_OPERATOR_TYPES \"")
        string(REPLACE ";" ";" OPERATOR_TYPES_STR "${OPERATOR_TYPES}")
        string(APPEND OUTPUT_CONTENT "${OPERATOR_TYPES_STR}\")\n\n")
    endif()
    
    # Template types
    if(TEMPLATE_TYPES)
        string(APPEND OUTPUT_CONTENT "set(DISCOVERED_TEMPLATE_TYPES \"")
        string(REPLACE ";" ";" TEMPLATE_TYPES_STR "${TEMPLATE_TYPES}")
        string(APPEND OUTPUT_CONTENT "${TEMPLATE_TYPES_STR}\")\n\n")
    endif()
    
    # All types (for reference)
    string(APPEND OUTPUT_CONTENT "set(ALL_DISCOVERED_TYPES \"")
    string(REPLACE ";" ";" ALL_TYPES_STR "${ALL_DISCOVERED_TYPES}")
    string(APPEND OUTPUT_CONTENT "${ALL_TYPES_STR}\")\n")
    
    # Write output file
    file(WRITE ${OUTPUT_FILE} "${OUTPUT_CONTENT}")
    
    # Report results
    list(LENGTH ALL_DISCOVERED_TYPES TOTAL_COUNT)
    list(LENGTH GLOBAL_TYPES GLOBAL_COUNT)
    list(LENGTH POPULATION_TYPES POP_COUNT)
    list(LENGTH OPERATOR_TYPES OP_COUNT)
    list(LENGTH TEMPLATE_TYPES TEMPLATE_COUNT)
    
    message(STATUS "Type discovery completed:")
    message(STATUS "  Total types discovered: ${TOTAL_COUNT}")
    message(STATUS "  Global scope types: ${GLOBAL_COUNT}")
    message(STATUS "  Population namespace types: ${POP_COUNT}")
    message(STATUS "  Operator namespace types: ${OP_COUNT}")
    message(STATUS "  Template types: ${TEMPLATE_COUNT}")
    message(STATUS "  Results saved to: ${OUTPUT_FILE}")
endfunction()