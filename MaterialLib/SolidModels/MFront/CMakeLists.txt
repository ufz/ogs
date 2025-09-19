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
    PRIVATE MathLib MeshLib ParameterLib
)
target_include_directories(
    MaterialLib_SolidModels_MFront PUBLIC ThirdParty/MGIS/include
)
target_compile_definitions(
    MaterialLib_SolidModels_MFront PRIVATE OGS_USE_MFRONT
)
if(NOT GUIX_BUILD AND NOT CONDA_BUILD)
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
    NeoHooke
    PowerLawLinearCreep
    SaintVenantKirchhoff
    StandardElasticityBrick
    StandardElasticityBrickOrtho
    ThermoPoroElasticity
    WeakPlaneModel
    WeakPlaneModel_nonIter
)

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
        ${CMAKE_COMMAND} -E env PATH="${TFELHOME}/bin:$ENV{PATH}"
                                CMAKE_POLICY_VERSION_MINIMUM="3.10"
                                --
        mfront
        --interface=generic
        --obuild=level0 ${_mfront_behaviour_files} # TODO: check for debug?
        --generator=cmake
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_generated_mfront_lib}
            ${_copied_mfront_lib}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND_EXPAND_LISTS
    DEPENDS ${_mfront_behaviour_files}
)
add_custom_target(build_mfront ALL DEPENDS ${_copied_mfront_lib})
add_dependencies(MaterialLib_SolidModels_MFront build_mfront)
install(FILES ${_copied_mfront_lib} TYPE LIB PERMISSIONS
    OWNER_READ GROUP_READ WORLD_READ
    OWNER_WRITE GROUP_WRITE
    OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
