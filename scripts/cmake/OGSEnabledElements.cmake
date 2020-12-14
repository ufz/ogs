# TODO: Move to process lib
# Definitions controlling which FEM elements will be compiled
if(NOT OGS_MAX_ELEMENT_DIM MATCHES "^[0-3]$")
  message(FATAL_ERROR "OGS_MAX_ELEMENT_DIM must be an integer between 0 and 3.")
endif()
add_definitions(-DOGS_MAX_ELEMENT_DIM=${OGS_MAX_ELEMENT_DIM})

if(NOT OGS_MAX_ELEMENT_ORDER MATCHES "^[0-9]$")
  message(FATAL_ERROR "OGS_MAX_ELEMENT_ORDER must be an integer.")
endif()
add_definitions(-DOGS_MAX_ELEMENT_ORDER=${OGS_MAX_ELEMENT_ORDER})


if(OGS_ENABLE_ELEMENT_SIMPLEX)
    add_definitions(-DOGS_ENABLE_ELEMENT_SIMPLEX)
endif()
if(OGS_ENABLE_ELEMENT_CUBOID)
    add_definitions(-DOGS_ENABLE_ELEMENT_CUBOID)
endif()
if(OGS_ENABLE_ELEMENT_PRISM)
    add_definitions(-DOGS_ENABLE_ELEMENT_PRISM)
endif()
if(OGS_ENABLE_ELEMENT_PYRAMID)
    add_definitions(-DOGS_ENABLE_ELEMENT_PYRAMID)
endif()
