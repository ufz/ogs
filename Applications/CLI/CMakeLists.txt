
add_executable(ogs
    ogs.cpp
)

target_link_libraries(ogs
    ApplicationsLib
    FileIO
)

if(OGS_BUILD_GUI)
	target_link_libraries(ogs
		QtDataView
		QtStratView
		VtkVis
		vtkInteractionStyle
		Qt4::QtCore
		Qt4::QtGui
		Qt4::QtXml
		Qt4::QtNetwork
	)
	if(CMAKE_CROSSCOMPILING)
		target_link_libraries(ogs
			${QT_XML_DEPS_LIBRARIES}
			${QT_GUI_DEPS_LIBRARIES}
			${QT_NETWORK_DEPS_LIBRARIES})
	endif()
endif()

ADD_CATALYST_DEPENDENCY(ogs)

if(OGS_USE_PETSC)
    target_link_libraries(ogs ${PETSC_LIBRARIES})
endif()

if(NOT OGS_USE_MPI)
    include(Tests.cmake)
endif()

####################
### Installation ###
####################
install(TARGETS ogs RUNTIME DESTINATION bin COMPONENT ogs_cli)

set(CPACK_PACKAGE_EXECUTABLES ${CPACK_PACKAGE_EXECUTABLES} "ogs" "OGS Simulator")
cpack_add_component(ogs_cli
	DISPLAY_NAME "OGS THMC Simulator"
	DESCRIPTION "The command line interface for OpenGeoSys."
	GROUP Applications
)

include(packaging/InstallDependencies)
InstallDependencies(ogs ogs_cli)
