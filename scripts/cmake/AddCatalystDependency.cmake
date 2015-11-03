macro(ADD_CATALYST_DEPENDENCY target)
	if(ParaView_FOUND)
		include("${PARAVIEW_USE_FILE}")

		# see http://stackoverflow.com/questions/18642155
		set_property(TARGET ${target} APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})
	elseif(VTK_FOUND)
		include( ${VTK_USE_FILE} )
	endif()

	if(TARGET VtkRescan)
		add_dependencies(${target} VtkRescan)
	endif()
endmacro()
