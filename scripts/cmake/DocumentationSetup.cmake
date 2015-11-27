if(DOXYGEN_FOUND)

	option(DOCS_GENERATE_DOCSET "Generate Dash Docsets." OFF)

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
		find_program(DOCSETUTIL_TOOLPATH docsetutil)
		if(NOT DOCSETUTIL_TOOLPATH)
			message(FATAL_ERROR "docsetutil required for Docset-generation!")
		endif()
		set(DOCS_GENERATE_TREEVIEW_STRING "NO" CACHE INTERNAL "")
		set(DOCS_DISABLE_INDEX_STRING "YES" CACHE INTERNAL "")
		set(DOCS_GENERATE_DOCSET_STRING "YES" CACHE INTERNAL "")
		set(DOCS_SEARCHENGINE_STRING "NO" CACHE INTERNAL "")
		add_custom_command(TARGET doc POST_BUILD
			COMMAND make WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/docs
			COMMENT "Generating docset ...")
	endif()

	configure_file(Documentation/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)

endif()
