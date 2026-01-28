# cmake-lint: disable=R0912,R0915,C0103

# Constants for header filtering patterns
set(_HEADER_EXCLUDE_PATTERNS
    ".*-impl\\.h$"
    ".*Dialog\\.h$"
    ".*Widget\\.h$"
    ".*Window\\.h$"
    "ui_.*\\.h$"
)

# Determine if a header file should be included in standalone compilation checks
function(_filter_header file_path out_var)
    # Only check header files
    if(NOT "${file_path}" MATCHES ".*\\.h$")
        set(${out_var} FALSE PARENT_SCOPE)
        return()
    endif()

    # Exclude implementation headers and Qt UI files
    foreach(pattern ${_HEADER_EXCLUDE_PATTERNS})
        if("${file_path}" MATCHES "${pattern}")
            set(${out_var} FALSE PARENT_SCOPE)
            return()
        endif()
    endforeach()

    set(${out_var} TRUE PARENT_SCOPE)
endfunction()

# Collect header files from a library target that should be checked
function(_collect_library_headers target out_var)
    # Early return if target doesn't exist
    if(NOT TARGET ${target})
        set(${out_var} "" PARENT_SCOPE)
        return()
    endif()

    get_target_property(SOURCE_FILES ${target} SOURCES)
    if(NOT SOURCE_FILES)
        set(${out_var} "" PARENT_SCOPE)
        return()
    endif()

    # Filter and collect headers
    set(headers_to_check "")
    foreach(file ${SOURCE_FILES})
        _filter_header("${file}" should_check)
        if(should_check)
            list(APPEND headers_to_check "${file}")
        endif()
    endforeach()

    set(${out_var} "${headers_to_check}" PARENT_SCOPE)
endfunction()

# Helper function to mangle header path into a valid target name
function(_mangle_header_name header_file out_name)
    string(REPLACE "${PROJECT_SOURCE_DIR}/" "" relative_path "${header_file}")
    string(REPLACE "/" "-" mangled_name "${relative_path}")
    string(REPLACE "." "_" mangled_name "${mangled_name}")
    set(${out_name} "${mangled_name}" PARENT_SCOPE)
endfunction()

# Create a compile-only check target for a single header
function(_create_header_check_target library_target header_file)
    # Mangle header path to create target name
    _mangle_header_name("${header_file}" mangled_name)
    set(target_name "check-header-${mangled_name}")

    # Create directory for generated sources and object files
    set(gen_dir "${CMAKE_BINARY_DIR}/CheckHeaderCompilation/${library_target}")
    file(MAKE_DIRECTORY "${gen_dir}")

    # Generate .cpp file
    set(gen_source "${gen_dir}/${mangled_name}.cpp")
    file(WRITE "${gen_source}"
        "// Auto-generated file to check standalone compilation\n"
        "#include \"${header_file}\"\n"
    )

    # Create the check target as an object library (compile-only, no linking)
    add_library(${target_name} OBJECT "${gen_source}")

    # Apply compilation settings
    target_include_directories(${target_name} PRIVATE
        $<TARGET_PROPERTY:${library_target},INTERFACE_INCLUDE_DIRECTORIES>
        ${PROJECT_SOURCE_DIR}
        $<$<BOOL:${OGS_BUILD_GUI}>:${PROJECT_SOURCE_DIR}/Applications/DataExplorer>
    )

    target_compile_definitions(${target_name} PRIVATE
        $<TARGET_PROPERTY:${library_target},INTERFACE_COMPILE_DEFINITIONS>
    )

    # Set target properties (C++ standard inherited globally)
    set_target_properties(${target_name} PROPERTIES
        EXCLUDE_FROM_ALL TRUE
        FOLDER "HeaderChecks"
    )

    # Add dependency on library target (for generated headers)
    if(TARGET ${library_target})
        add_dependencies(${target_name} ${library_target})
    endif()
endfunction()

# Create check targets for all headers in a library
function(check_header_compilation_for_library library_target)
    # Early return if target doesn't exist
    if(NOT TARGET ${library_target})
        message(STATUS "Target ${library_target} does not exist, skipping header checks.")
        return()
    endif()

    # Collect headers to check
    _collect_library_headers(${library_target} headers_to_check)

    if(NOT headers_to_check)
        message(STATUS "No headers to check for ${library_target}")
        return()
    endif()

    # Get source directory for status message
    get_target_property(SOURCE_DIR ${library_target} SOURCE_DIR)
    string(REPLACE "${PROJECT_SOURCE_DIR}/" "" DIRECTORY ${SOURCE_DIR})

    list(LENGTH headers_to_check num_headers)
    message(STATUS "Setting up header compilation checks for ${DIRECTORY} (${num_headers} headers)...")

    # Create check target for each header
    set(all_header_targets "")
    foreach(header ${headers_to_check})
        _create_header_check_target(${library_target} "${header}")

        # Calculate target name using the same mangling logic
        _mangle_header_name("${header}" mangled_name)
        set(target_name "check-header-${mangled_name}")

        list(APPEND all_header_targets ${target_name})
    endforeach()

    # Create library-level phony target
    set(library_check_target "check-headers-${library_target}")
    add_custom_target(${library_check_target})

    if(all_header_targets)
        add_dependencies(${library_check_target} ${all_header_targets})
    endif()
endfunction()

# Core libraries to check for header compilation
set(_CORE_LIBS_TO_CHECK
    BaseLib
    ChemistryLib
    GeoLib
    GitInfoLib
    CMakeInfoLib
    TestInfoLib
    MaterialLib
    MathLib
    MeshGeoToolsLib
    MeshLib
    NumLib
    ParameterLib
    ProcessLib
    ApplicationsLib
    ApplicationsFileIO
    DataHolderLib
)

# GUI libraries to check (only when GUI is enabled)
set(_GUI_LIBS_TO_CHECK
    QtBase
    QtDataView
    QtDiagramView
    VtkVis
    QtStratView
)

# Main function: Create header compilation check targets for all libraries
function(check_header_compilation)
    if(NOT OGS_CHECK_HEADER_COMPILATION)
        return()
    endif()

    message(STATUS "")
    message(STATUS "========================================")
    message(STATUS "Setting up header compilation checks")
    message(STATUS "========================================")

    # Build list of libraries to check
    set(LIBS_TO_CHECK ${_CORE_LIBS_TO_CHECK})
    if(OGS_BUILD_GUI)
        list(APPEND LIBS_TO_CHECK ${_GUI_LIBS_TO_CHECK})
    endif()

    # Create check targets for each library
    set(all_library_targets "")
    foreach(lib ${LIBS_TO_CHECK})
        check_header_compilation_for_library(${lib})
        set(lib_target "check-headers-${lib}")
        if(TARGET ${lib_target})
            list(APPEND all_library_targets ${lib_target})
        endif()
    endforeach()

    # Create top-level target
    add_custom_target(check-headers)
    if(all_library_targets)
        add_dependencies(check-headers ${all_library_targets})
    endif()

    message(STATUS "")
    message(STATUS "Header compilation check targets created.")
    message(STATUS "Run 'ninja check-headers' to execute all checks.")
    message(STATUS "Run 'ninja check-headers-<library>' to check a specific library.")
    message(STATUS "========================================")
    message(STATUS "")
endfunction()
