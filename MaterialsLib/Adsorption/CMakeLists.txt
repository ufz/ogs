file(GLOB Adsorption_SOURCES *.cpp)
file(GLOB Adsorption_HEADERS *.h)

add_library(MaterialsLibAdsorption ${Adsorption_HEADERS} ${Adsorption_SOURCES} )
