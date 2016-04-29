file(GLOB Adsorption_SOURCES *.cpp)
file(GLOB Adsorption_HEADERS *.h)

add_library(MaterialsLibAdsorption STATIC ${Adsorption_HEADERS} ${Adsorption_SOURCES} )
