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

target_link_libraries(
    MaterialLib_SolidModels_MFront PUBLIC BaseLib NumLib MFrontGenericInterface
                                          Boost::boost PRIVATE MathLib MeshLib
                                                               ParameterLib
)
target_include_directories(
    MaterialLib_SolidModels_MFront PUBLIC ThirdParty/MGIS/include
)
target_compile_definitions(
    MaterialLib_SolidModels_MFront PRIVATE OGS_USE_MFRONT
)
if(NOT GUIX_BUILD)
    install(TARGETS MFrontGenericInterface)
endif()

# Setup / build mfront behaviours lib
set(_mfront_behaviours
    BDT
    DruckerPrager
    GuentherSalzer
    GuentherSalzerDilatancy_semi_expl
    HoekBrownC2
    Lubby2
    Lubby2mod
    ModCamClay_semiExpl
    ModCamClay_semiExpl_absP
    ModCamClay_semiExpl_constE
    MohrCoulombAbboSloan
    MohrCoulombAbboSloanAniso
    MohrCoulombAbboSloanOrtho
    MohrCoulombAbboSloanUBI
    MohrCoulombAbboSloanUBIOrtho
    PowerLawLinearCreep
    SaintVenantKirchhoff
    StandardElasticityBrick
    StandardElasticityBrickOrtho
    ThermoPoroElasticity
)

if((OGS_BUILD_WHEEL AND NOT APPLE_ARM) OR GUIX_BUILD)
    # Build behaviours lib with mfront executable
    foreach(behaviour ${_mfront_behaviours})
        list(APPEND _mfront_behaviour_files
             "${CMAKE_CURRENT_SOURCE_DIR}/${behaviour}.mfront"
        )
    endforeach()
    set(_generated_mfront_lib
        ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
    set(_copied_mfront_lib
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libOgsMFrontBehaviour${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
    add_custom_command(
        OUTPUT ${_copied_mfront_lib}
        COMMAND
            ${CMAKE_COMMAND} -E env PATH="${TFELHOME}/bin:$ENV{PATH}" mfront
            --interface=generic --silent-build=true # suppresses output
                                                    # (compiler warnings), does
                                                    # not work?
            --obuild=level0 ${_mfront_behaviour_files} # TODO: check for debug?
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_generated_mfront_lib}
                ${_copied_mfront_lib}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND_EXPAND_LISTS
        DEPENDS ${_mfront_behaviour_files}
    )
    add_custom_target(
        InternalOgsMFrontBehaviour ALL DEPENDS ${_copied_mfront_lib}
    )
    add_library(OgsMFrontBehaviour SHARED IMPORTED)
    set_target_properties(
        OgsMFrontBehaviour PROPERTIES IMPORTED_LOCATION ${_copied_mfront_lib}
                                      IMPORTED_NO_SONAME TRUE
    )
    install(FILES ${_copied_mfront_lib} TYPE LIB)
else()
    # Build behaviours lib as part of ogs' build process
    mfront_behaviours_check_library(OgsMFrontBehaviour ${_mfront_behaviours})
    # ------------------------------ ATTENTION --------------------------------
    # Disable warnings for generated OgsMFrontBehaviour. To enable warnings
    # comment out the following target_compile_options()-statement:
    target_compile_options(
        OgsMFrontBehaviour PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                                   $<$<CXX_COMPILER_ID:MSVC>:/W0>
    )
    target_link_libraries(
        MaterialLib_SolidModels_MFront PUBLIC OgsMFrontBehaviour
    )
    install(TARGETS OgsMFrontBehaviour)
endif()
