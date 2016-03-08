file(GLOB adsorption_SOURCES *.cpp)
file(GLOB adsorption_HEADERS *.h)

# set( HEADERS
# 	adsorption.h
# 	density_legacy.h
# 	density_hauer.h
# 	density_mette.h
# )
#
# set( SOURCES
# 	adsorption.cpp
# 	density_legacy.cpp
# 	density_hauer.cpp
# 	density_mette.cpp
# )


# include_directories(
# 	${CMAKE_SOURCE_DIR}/FEM
# )

add_library(MaterialsLib_adsorption STATIC ${adsorption_HEADERS} ${adsorption_SOURCES} )
