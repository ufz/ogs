cmake_policy(SET CMP0009 NEW)

# Moves an *.md file from the src directory to DocAux directory in the build
# tree augmenting it with basic extra information, such as the list of child
# pages and the page title.
function(documentationProjectFilePutIntoPlace p)
    file(RELATIVE_PATH relative_path ${DocumentationProjectFileInputDir} ${p})
    get_filename_component(dir_name ${relative_path} DIRECTORY)

    get_filename_component(otagname ${relative_path} NAME_WE)
    if (otagname MATCHES ^[ic]_)
        # if the file name starts with an underscore, then this files is
        # the "table of contents of the current directory

        file(MAKE_DIRECTORY "${DocumentationProjectFileBuildDir}/${dir_name}")

        set(postfix "# Child parameters, attributes and cases\n\n")

        # gather other parameter files
        # the loop below will effects a page hierarchy to be built
        file(GLOB param_files ${DocumentationProjectFileInputDir}/${dir_name}/*)
        set(subpagelist "")
        foreach(pf ${param_files})
            # ignore hidden files
            if (pf MATCHES /[.][^/]+)
                continue()
            endif()

            get_filename_component(rel_pf ${pf} NAME_WE)

            # if the file name matches ^[ic]_, then this
            # is the "table of contents" file already processed outside
            # of this loop
            if (NOT rel_pf MATCHES ^[ic]_)
                if(IS_DIRECTORY ${pf})
                    set(pf_tagname ${rel_pf})
                else()
                    if (NOT "${rel_pf}" MATCHES ^._)
                        message(SEND_ERROR "Path `${rel_pf}' has a wrong name."
                            " Full path is `${pf}'.")
                        continue()
                    endif()

                    string(SUBSTRING "${rel_pf}" 2 -1 pf_tagname)
                endif()

                if ("${dir_name}" STREQUAL "") # toplevel dir must be treated slightly different
                    set(pf_tagpath "${pf_tagname}")
                else()
                    set(pf_tagpath "${dir_name}/${pf_tagname}")
                    string(REPLACE "/" "__" pf_tagpath "${pf_tagpath}")
                endif()
                message("  t.o.c. entry ${pf_tagpath}")

                if (rel_pf MATCHES ^a_)
                    set(pagenameprefix "ogs_file_attr__")
                else()
                    set(pagenameprefix "ogs_file_param__")
                endif()

                list(FIND subpagelist "${pagenameprefix}${pf_tagpath}" idx)
                if (NOT idx EQUAL -1)
                    message(SEND_ERROR "The subpagelist already contains"
                        " ${pagenameprefix}${pf_tagpath}. Maybe there are"
                        " duplicate documentation files.")
                else()
                    list(APPEND subpagelist "${pagenameprefix}${pf_tagpath}")
                endif()

                if (NOT IS_DIRECTORY "${pf}")
                    documentationProjectFilePutIntoPlace("${pf}")
                endif()
            endif()
        endforeach()

        list(SORT subpagelist)
        foreach(subpage ${subpagelist})
            set(postfix "${postfix} - \\subpage ${subpage}\n")
        endforeach()
    else()
        set(postfix "")
    endif()

    string(SUBSTRING ${otagname} 2 -1 tagname)
    if (dir_name STREQUAL "") # toplevel dir must be treated slightly different
        set(tagpath "${tagname}")
    else()
        if (otagname MATCHES ^[ic]_) # treat "table of contents" file special
            string(REPLACE "/" "__" tagpath "${dir_name}")
        else()
            string(REPLACE "/" "__" tagpath "${dir_name}/${tagname}")
        endif()
    endif()
    message("  child param  ${tagpath}")

    set(pagenameprefix "ogs_file_param__")
    if (otagname MATCHES ^i_ AND dir_name STREQUAL "")
        set(pagetitle "OGS Input File Parameters")
    elseif(otagname MATCHES ^c_)
        set(pagetitle "[case]&emsp;${tagname}")
    elseif(otagname MATCHES ^t_ OR otagname MATCHES ^i_)
        set(pagetitle "[tag]&emsp;${tagname}")
    elseif(otagname MATCHES ^a_)
        set(pagetitle "[attr]&emsp;${tagname}")
        set(pagenameprefix "ogs_file_attr__")
    else()
        message(SEND_ERROR "Tag name ${otagname} does not match in any case."
            " Maybe there is a file with a wrong name in the documentation"
            " directory.")
    endif()

    # read, augment, write file content
    file(READ ${p} content)
    set(content "/*! \\page ${pagenameprefix}${tagpath} ${pagetitle}\n${content}\n\n${postfix}\n")
    if (NOT doc_use_external_tools)
        set(ending "\n*/\n")
    else()
        set(ending "") # external tools shall finish the file
    endif()
    string(REGEX REPLACE .md$ .dox output_file "${DocumentationProjectFileBuildDir}/${relative_path}")
    file(WRITE "${output_file}" "${content}${ending}")
endfunction()


set(DocumentationProjectFileBuildDir ${PROJECT_BINARY_DIR}/DocAux/dox/ProjectFile)
set(DocumentationProjectFileInputDir ${PROJECT_SOURCE_DIR}/Documentation/ProjectFile)

# remove old output
if (IS_DIRECTORY ${DocumentationProjectFileBuildDir})
    file(REMOVE_RECURSE ${DocumentationProjectFileBuildDir})
endif()

# traverse input file hierarchy
file(GLOB_RECURSE input_paths FOLLOW_SYMLINKS
    ${DocumentationProjectFileInputDir}/c_* ${DocumentationProjectFileInputDir}/i_*)

foreach(p ${input_paths})
    message("directory index file ${p}")
    documentationProjectFilePutIntoPlace(${p})
endforeach()
