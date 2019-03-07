set(SOURCES CreateMFront.cpp CreateMFront.h)

if(OGS_USE_MFRONT)
    list(APPEND SOURCES MFront.cpp MFront.h)
endif()

add_library(MaterialLib_SolidModels_MFront ${SOURCES})

if(BUILD_SHARED_LIBS)
    install(TARGETS MaterialLib_SolidModels_MFront LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                                   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

target_link_libraries(MaterialLib_SolidModels_MFront
    PUBLIC BaseLib NumLib logog
    PRIVATE MathLib MeshLib
)

if (OGS_USE_MFRONT)
    target_include_directories(MaterialLib_SolidModels_MFront
        PUBLIC ${MGIS_INCLUDE_DIR}
    )
    target_link_libraries(MaterialLib_SolidModels_MFront
        PUBLIC ${MGIS_LIBRARY}
    )
    target_compile_definitions(MaterialLib_SolidModels_MFront PRIVATE OGS_USE_MFRONT)
endif()
