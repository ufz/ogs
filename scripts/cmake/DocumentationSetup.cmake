if(DOXYGEN_FOUND)

    option(DOCS_GENERATE_DOCSET "Generate Dash Docsets." OFF)
    option(DOCS_GENERATE_LOGFILE "Outputs Doxygen warnings to a file in the build directory." OFF)

    set(DOT_FOUND "NO")
    if(DOXYGEN_DOT_FOUND)
        set(DOT_FOUND "YES")
    endif()

    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating source code documentation with Doxygen." VERBATIM)

    # Defaults
    set(DOCS_GENERATE_TREEVIEW_STRING "YES" CACHE INTERNAL "")
    set(DOCS_DISABLE_INDEX_STRING "NO" CACHE INTERNAL "")
    set(DOCS_GENERATE_DOCSET_STRING "NO" CACHE INTERNAL "")
    set(DOCS_SEARCHENGINE_STRING "YES" CACHE INTERNAL "")

    # Dash Docsets
    if(DOCS_GENERATE_DOCSET)
        find_program(DOCSETUTIL_TOOLPATH docsetutil
            PATH /Applications/Xcode.app/Contents/Developer/usr/bin)
        if(NOT DOCSETUTIL_TOOLPATH)
            message(FATAL_ERROR "docsetutil required for Docset-generation!")
        endif()
        set(DOCS_GENERATE_TREEVIEW_STRING "NO" CACHE INTERNAL "")
        set(DOCS_DISABLE_INDEX_STRING "YES" CACHE INTERNAL "")
        set(DOCS_GENERATE_DOCSET_STRING "YES" CACHE INTERNAL "")
        set(DOCS_SEARCHENGINE_STRING "NO" CACHE INTERNAL "")
        add_custom_command(TARGET doc POST_BUILD
            COMMAND make
            COMMAND mv org.doxygen.Project.docset ogs6.docset
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/docs
            COMMENT "Generating docset ...")
        configure_file(Documentation/DocsetFeed.xml.in ${PROJECT_BINARY_DIR}/docs/ogs6.xml)
    endif()

    if(DOCS_GENERATE_LOGFILE)
        set(OGS_DOXYGEN_LOGFILE "${PROJECT_BINARY_DIR}/DoxygenWarnings.log" CACHE INTERNAL "")
    endif()

    set(DOXYGEN_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/docs CACHE INTERNAL "")

    if(OGS_WEB_EMBED_DOXYGEN)
        set(DOXYGEN_OUTPUT_DIRECTORY
            ${PROJECT_SOURCE_DIR}/web/public/docs/doxygen
            CACHE INTERNAL "")
        set(DOXYGEN_HTML_HEADER
            ${PROJECT_SOURCE_DIR}/web/public/docs/doxygen/header/index.html
            CACHE INTERNAL "")
        set(DOXYGEN_HTML_FOOTER
            ${PROJECT_SOURCE_DIR}/web/public/docs/doxygen/footer/index.html
            CACHE INTERNAL "")
        set(DOXYGEN_HTML_EXTRA_STYLESHEET
            ${PROJECT_SOURCE_DIR}/web/static/css/doxygen-extra.css)
        add_dependencies(doc web)

        if(NPM AND GRUNT)
            add_custom_command(TARGET doc POST_BUILD
                COMMAND ${GRUNT}
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/web
                COMMENT "Running grunt for prefix-ing doxygen css...")
        endif()
    endif()

    configure_file(Documentation/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)

    if (BASH_TOOL_PATH AND PYTHON_EXECUTABLE)
        set(doc_use_external_tools TRUE)
    else()
        set(doc_use_external_tools FALSE)
    endif()

    # TODO that will always transform all of the input files no matter if they changed
    # maybe this behaviour can be changed to on-demand processing
    add_custom_target(internal_pre_doc
        ${CMAKE_COMMAND}
        -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
        -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        -Ddoc_use_external_tools=${doc_use_external_tools}
        -P ${PROJECT_SOURCE_DIR}/scripts/cmake/DocumentationProjectFile.cmake
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating project file documentation hierarchy." VERBATIM)
    add_dependencies(doc internal_pre_doc)

    if (doc_use_external_tools)
        set(data_dir "${Data_SOURCE_DIR}")
        add_custom_target(internal_pre_doc_qa_page
            ${BASH_TOOL_PATH}
            "${PROJECT_SOURCE_DIR}/scripts/doc/generate-project-file-doc-qa.sh"
            ${PROJECT_SOURCE_DIR}
            ${PROJECT_BINARY_DIR}
            ${data_dir}
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Generating project file documentation quality assurance pages." VERBATIM)
        add_dependencies(doc internal_pre_doc_qa_page)
        add_dependencies(internal_pre_doc_qa_page internal_pre_doc)
    endif()
endif()
