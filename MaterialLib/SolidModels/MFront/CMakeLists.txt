set(SOURCES
    CreateMFront.cpp
    CreateMFront.h
    MFront.cpp
    MFront.h
    MFrontGeneric.cpp
    MFrontGeneric.h
    CreateMFrontGeneric.cpp
    CreateMFrontGeneric.h
)

ogs_add_library(MaterialLib_SolidModels_MFront ${SOURCES})

mfront_behaviours_check_library(
    OgsMFrontBehaviour
    BDT
    DruckerPrager
    GuentherSalzer
    Lubby2
    Lubby2mod
    ModCamClay_semiExpl_constE
    ModCamClay_semiExplParaInitNLnu_abs
    ModCamClay_semiExplParaInitNLnu_inc
    MohrCoulombAbboSloan
    MohrCoulombAbboSloanAniso
    MohrCoulombAbboSloanOrtho
    MohrCoulombAbboSloanUBI
    MohrCoulombAbboSloanUBIOrtho
    PowerLawLinearCreep
    StandardElasticityBrick
    StandardElasticityBrickOrtho
    ThermoPoroElasticity
)

target_link_libraries(
    MaterialLib_SolidModels_MFront PUBLIC BaseLib NumLib MFrontGenericInterface
                                          OgsMFrontBehaviour Boost::boost
    PRIVATE MathLib MeshLib ParameterLib
)

# Disable warnings for generated OgsMFrontBehaviour
target_compile_options(
    OgsMFrontBehaviour PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                               $<$<CXX_COMPILER_ID:MSVC>:/W0>
)

target_include_directories(
    MaterialLib_SolidModels_MFront PUBLIC ThirdParty/MGIS/include
)
target_compile_definitions(
    MaterialLib_SolidModels_MFront PRIVATE OGS_USE_MFRONT
)

install(TARGETS MFrontGenericInterface OgsMFrontBehaviour)
