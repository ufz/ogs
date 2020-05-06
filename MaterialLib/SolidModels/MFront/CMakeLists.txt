set(SOURCES CreateMFront.cpp CreateMFront.h)

list(APPEND SOURCES MFront.cpp MFront.h)

ogs_add_library(MaterialLib_SolidModels_MFront ${SOURCES})

target_link_libraries(MaterialLib_SolidModels_MFront
                      PUBLIC BaseLib NumLib OgsMFrontBehaviour
                      PRIVATE MathLib MeshLib spdlog::spdlog)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/ThirdParty/MGIS")
include(cmake/modules/tfel)
mfront_behaviours_check_library(OgsMFrontBehaviour
                                BDT
                                DruckerPrager
                                Elasticity
                                MohrCoulombAbboSloan
                                MohrCoulombAbboSloanOrtho
                                StandardElasticityBrick
                                StandardElasticityBrickOrtho
                                Lubby2
                                Lubby2mod
                                MohrCoulombAbboSloanAniso
                                MohrCoulombAbboSloanUBI
                                MohrCoulombAbboSloanUBIOrtho)

target_include_directories(MaterialLib_SolidModels_MFront
                           PUBLIC ThirdParty/MGIS/include)
target_link_libraries(MaterialLib_SolidModels_MFront
                      PUBLIC MFrontGenericInterface)
target_compile_definitions(MaterialLib_SolidModels_MFront PRIVATE
                           OGS_USE_MFRONT)
