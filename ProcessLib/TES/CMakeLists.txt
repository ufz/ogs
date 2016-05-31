if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl-dynamic.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl.h" COPYONLY
        )
else()
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl-fixed.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/TESLocalAssemblerInner-impl-incl.h" COPYONLY
        )
endif()
