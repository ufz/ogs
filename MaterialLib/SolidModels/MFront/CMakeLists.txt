set(SOURCES CreateMFront.cpp CreateMFront.h)

list(APPEND SOURCES MFront.cpp MFront.h)

add_library(MaterialLib_SolidModels_MFront ${SOURCES})

if(BUILD_SHARED_LIBS)
    install(TARGETS MaterialLib_SolidModels_MFront
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

target_link_libraries(MaterialLib_SolidModels_MFront
                      PUBLIC BaseLib NumLib logog OgsMFrontBehaviour
                      PRIVATE MathLib MeshLib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/ThirdParty/MGIS")
include(cmake/modules/tfel)
mfront_behaviours_check_library(OgsMFrontBehaviour Elasticity
                                StandardElasticityBrick)

target_include_directories(MaterialLib_SolidModels_MFront
                           PUBLIC ThirdParty/MGIS/include)
target_link_libraries(MaterialLib_SolidModels_MFront
                      PUBLIC MFrontGenericInterface)
target_compile_definitions(MaterialLib_SolidModels_MFront PRIVATE
                           OGS_USE_MFRONT)
