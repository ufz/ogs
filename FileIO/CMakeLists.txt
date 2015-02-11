set(SOURCES
	AsciiRasterInterface.h
	AsciiRasterInterface.cpp
	GMSInterface.h
	GMSInterface.cpp
	GMSHInterface.h
	GMSHInterface.cpp
	PetrelInterface.h
	PetrelInterface.cpp
	readMeshFromFile.h
	readMeshFromFile.cpp
	SHPInterface.h
	SHPInterface.cpp
	TetGenInterface.h
	TetGenInterface.cpp
	TINInterface.h
	TINInterface.cpp
	Writer.h
	Writer.cpp
	writeMeshToFile.h
	writeMeshToFile.cpp
)

GET_SOURCE_FILES(SOURCES_LEGACY Legacy)
GET_SOURCE_FILES(SOURCES_GMSHIO GmshIO)
GET_SOURCE_FILES(SOURCES_RAPID_XML RapidXmlIO)
GET_SOURCE_FILES(SOURCES_VTK VtkIO)
set(SOURCES ${SOURCES} ${SOURCES_LEGACY} ${SOURCES_GMSHIO} ${SOURCES_RAPID_XML} ${SOURCES_VTK})

GET_SOURCE_FILES(SOURCES_BASE_XML XmlIO)
set(SOURCES ${SOURCES} ${SOURCES_BASE_XML})
GET_SOURCE_FILES(SOURCES_BOOST_XML XmlIO/Boost)
set(SOURCES ${SOURCES} ${SOURCES_BOOST_XML})
if(QT4_FOUND)
	set(SOURCES ${SOURCES} FEFLOWInterface.h FEFLOWInterface.cpp)
	GET_SOURCE_FILES(SOURCES_QT_XML XmlIO/Qt)
	set(SOURCES ${SOURCES} ${SOURCES_QT_XML})
endif()

# It could be used for other MPI based DDC approach in future.
if(OGS_USE_PETSC)
	GET_SOURCE_FILES(SOURCES_MPI_IO MPI_IO)
	set(SOURCES ${SOURCES} ${SOURCES_MPI_IO})
endif()

# Create the library
add_library(FileIO STATIC ${SOURCES})

include_directories(
	.
	${CMAKE_CURRENT_SOURCE_DIR}/../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib
	${CMAKE_CURRENT_BINARY_DIR}/../BaseLib
)

target_link_libraries(FileIO
	GeoLib
	MeshLib
	zlib
	shp
)

ADD_CATALYST_DEPENDENCY(FileIO)

file(GLOB XSD_FILES . *.xsd)
if(APPLE AND OGS_BUILD_GUI)
	install(FILES ${XSD_FILES} DESTINATION ${CMAKE_BINARY_DIR}/_CPack_Packages/Darwin/DragNDrop/${CPACK_PACKAGE_FILE_NAME}/ALL_IN_ONE/DataExplorer.app/Contents/MacOS COMPONENT ogs_gui)
else()
	install(FILES ${XSD_FILES} DESTINATION bin COMPONENT ogs_cli)
endif()
