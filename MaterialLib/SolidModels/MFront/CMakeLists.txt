# Setup / build mfront behaviours lib
set(_mfront_behaviours
    BDT
    DruckerPrager
    GuentherSalzer
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
    StandardElasticityBrick
    StandardElasticityBrickOrtho
    ThermoPoroElasticity
)

if(NOT OGS_BUILD_WHEEL)
    mfront_behaviours_check_library(OgsMFrontBehaviour ${_mfront_behaviours})
    # Disable warnings for generated OgsMFrontBehaviour
    target_compile_options(
        OgsMFrontBehaviour PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-w>
                                   $<$<CXX_COMPILER_ID:MSVC>:/W0>
    )
else()
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
            --obuild ${_mfront_behaviour_files} # TODO: check for debug?
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
endif()

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
    MaterialLib_SolidModels_MFront
    PUBLIC BaseLib NumLib MFrontGenericInterface
           $<$<NOT:$<BOOL:${OGS_BUILD_WHEEL}>>:OgsMFrontBehaviour> Boost::boost
    PRIVATE MathLib MeshLib ParameterLib
)

target_include_directories(
    MaterialLib_SolidModels_MFront PUBLIC ThirdParty/MGIS/include
)
target_compile_definitions(
    MaterialLib_SolidModels_MFront PRIVATE OGS_USE_MFRONT
)

install(TARGETS MFrontGenericInterface)
