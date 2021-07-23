if(NOT DOXYGEN_FOUND)
    return()
endif()

# Fix for https://github.com/doxygen/doxygen/issues/6725
set(DOXYGEN_LATEX_MAKEINDEX_CMD "makeindex")

set(DOXYGEN_EXCLUDE
    ${PROJECT_SOURCE_DIR}/ThirdParty ${PROJECT_SOURCE_DIR}/scripts
    ${PROJECT_SOURCE_DIR}/Tests ${PROJECT_SOURCE_DIR}/Documentation/ProjectFile
    ${PROJECT_SOURCE_DIR}/web
)
set(DOXYGEN_FILE_PATTERNS *.h *.cpp *.tpp *.dox)
set(DOXYGEN_PROJECT_NAME "OGS")
set(DOXYGEN_PROJECT_NUMBER "${OGS_GIT_BRANCH}")
set(DOXYGEN_PROJECT_VERSION ${GIT_SHA1_SHORT})
set(DOXYGEN_PROJECT_LOGO
    ${PROJECT_SOURCE_DIR}/Documentation/OpenGeoSys-Logo.png
)
set(DOXYGEN_HTML_OUTPUT ${PROJECT_BINARY_DIR}/docs)
if("${OGS_GIT_BRANCH}" MATCHES "^v[0-9]\\.[0-9]\\.[0-9]")
    set(DOXYGEN_HTML_COLORSTYLE_HUE 190)
endif()
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_EXTRACT_PRIVATE YES)
set(DOXYGEN_EXTRACT_PRIV_VIRTUAL YES)
set(DOXYGEN_EXTRACT_PACKAGE YES)
set(DOXYGEN_EXTRACT_STATIC YES)
set(DOXYGEN_EXTRACT_ANON_NSPACES YES)
set(DOXYGEN_INTERNAL_DOCS YES)
set(DOXYGEN_SORT_BY_SCOPE_NAME YES)
set(DOXYGEN_SORT_MEMBERS_CTORS_1ST YES)
set(DOXYGEN_LAYOUT_FILE ${PROJECT_SOURCE_DIR}/Documentation/DoxygenLayout.xml)
set(DOXYGEN_CITE_BIB_FILES
    ${PROJECT_SOURCE_DIR}/Documentation/bibliography/ogs
    ${PROJECT_SOURCE_DIR}/Documentation/bibliography/other
)
set(DOXYGEN_QUIET YES)
set(DOXYGEN_WARN_LOGFILE ${PROJECT_BINARY_DIR}/DoxygenWarnings.log)
set(DOXYGEN_EXCLUDE_PATTERNS moc_* ui_* CMake*)
set(DOXYGEN_IMAGE_PATH ${PROJECT_SOURCE_DIR}
                       ${PROJECT_SOURCE_DIR}/Documentation/images
)
set(DOXYGEN_SOURCE_BROWSER YES)
set(DOXYGEN_INLINE_SOURCES YES)
set(DOXYGEN_REFERENCED_BY_RELATION YES)
set(DOXYGEN_REFERENCES_RELATION YES)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_MATHJAX_RELPATH
    https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/
)
set(DOXYGEN_GENERATE_LATEX NO)
set(DOXYGEN_EXTRA_PACKAGES amsmath amsfonts)
set(DOXYGEN_PREDEFINED
    CVODE_FOUND
    DOXYGEN_DOCU_ONLY
    GEOTIFF_FOUND
    OGS_BUILD_GUI
    OGS_ENABLE_ELEMENT_CUBOID
    OGS_ENABLE_ELEMENT_PRISM
    OGS_ENABLE_ELEMENT_PYRAMID
    OGS_ENABLE_ELEMENT_SIMPLEX
    OGS_USE_MFRONT
    OGS_USE_NETCDF
    OGS_USE_PYTHON
    USE_INSITU
    USE_LIS
    USE_MKL
    USE_PETSC
)

set(DOXYGEN_DOT_IMAGE_FORMAT "svg")
set(DOXYGEN_INTERACTIVE_SVG YES)

set(DOXYGEN_ALIASES
    "per{1} = \\1<sup>-1</sup>"
    "ogs_file_param{1} = \\xrefitem ogs_file_param \\\"Input File Parameter\\\" \\\"List of all Input File Parameters\\\" \\ref ogs_file_param__\\1 \\\"\\1\\\""
    "ogs_file_attr{1} = \\xrefitem ogs_file_param \\\"Input File Parameter\\\" \\\"List of all Input File Parameters\\\" \\ref ogs_file_attr__\\1 \\\"\\1\\\""
    "ogs_file_special = \\xrefitem ogs_file_param \\\"Input File Parameter\\\" \\\"List of all Input File Parameters\\\" special OGS input file parameter"
    "ogs_file_param_special{1} = \\xrefitem ogs_file_param \\\"Input File Parameter\\\" \\\"List of all Input File Parameters\\\" \\ref ogs_file_param__\\1 \\\"\\1\\\""
    "ogs_file_attr_special{1} = \\xrefitem ogs_file_param \\\"Input File Parameter\\\" \\\"List of all Input File Parameters\\\" \\ref ogs_file_attr__\\1 \\\"\\1\\\""
    "ogs_missing_documentation = \\xrefitem ogs_missing_documentation \\\"Missing Documentation\\\" \\\"OGS Input File Parameters&mdash\;List of incomplete documentation pages\\\" Documentation missing/incomplete"
)
set(DOXYGEN_CREATE_SUBDIRS YES)
if($ENV{CI_COMMIT_BRANCH} MATCHES "master|^v[0-9]\.[0-9]\.[0-9]")
    set(DOXYGEN_SERVER_BASED_SEARCH YES)
    set(DOXYGEN_EXTERNAL_SEARCH YES)
    set(DOXYGEN_SEARCHENGINE_URL
        "https://doxysearch.opengeosys.org/$ENV{CI_COMMIT_BRANCH}/doxysearch.cgi"
    )
    message(STATUS "Doxygen search server: ${DOXYGEN_SEARCHENGINE_URL}")
endif()
if(EXISTS ${PROJECT_BINARY_DIR}/cpp-dependencies.svg)
    set(_subpage_dependencies "* \\subpage Dependencies")
endif()
configure_file(
    ${PROJECT_SOURCE_DIR}/Documentation/mainpage.dox.in
    ${PROJECT_BINARY_DIR}/DocAux/dox/mainpage.dox
)

doxygen_add_docs(doc ${PROJECT_SOURCE_DIR}/ ${PROJECT_BINARY_DIR}/DocAux/dox)

if(BASH_TOOL_PATH AND Python3_EXECUTABLE)
    set(_doc_use_external_tools TRUE)
else()
    set(_doc_use_external_tools FALSE)
endif()
# TODO that will always transform all of the input files no matter if they
# changed maybe this behaviour can be changed to on-demand processing
add_custom_target(
    internal_pre_doc
    ${CMAKE_COMMAND}
    -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
    -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -Ddoc_use_external_tools=${_doc_use_external_tools}
    -P
    ${PROJECT_SOURCE_DIR}/scripts/cmake/DocumentationProjectFile.cmake
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating project file documentation hierarchy."
    VERBATIM
)
add_dependencies(doc internal_pre_doc)

if(_doc_use_external_tools)
    add_custom_target(
        internal_pre_doc_qa_page
        ${BASH_TOOL_PATH}
        "${PROJECT_SOURCE_DIR}/scripts/doc/generate-project-file-doc-qa.sh"
        ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} "${Data_SOURCE_DIR}"
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating project file documentation quality assurance pages."
        VERBATIM
    )
    add_dependencies(doc internal_pre_doc_qa_page)
    add_dependencies(internal_pre_doc_qa_page internal_pre_doc)
endif()
