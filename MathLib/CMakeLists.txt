# Source files grouped by a directory
get_source_files(SOURCES)

append_source_files(SOURCES Curve)
append_source_files(SOURCES InterpolationAlgorithms)
append_source_files(SOURCES Integration)
append_source_files(SOURCES LinAlg)
append_source_files(SOURCES LinAlg/Dense)
if(CVODE_FOUND)
    append_source_files(SOURCES ODE)
endif()
append_source_files(SOURCES Nonlinear)

if(OGS_USE_EIGEN)
    append_source_files(SOURCES LinAlg/Eigen)
endif()

if(OGS_USE_LIS)
    append_source_files(SOURCES LinAlg/Lis)
    append_source_files(SOURCES LinAlg/EigenLis)
endif()

if(OGS_USE_PETSC)
    append_source_files(SOURCES LinAlg/PETSc)
endif()

# Create the library
ogs_add_library(MathLib ${SOURCES})
target_include_directories(MathLib
    PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        $<$<BOOL:${OGS_USE_PETSC}>:${PETSC_INCLUDE_DIRS}>
    PRIVATE
        $<$<BOOL:${OGS_USE_CVODE}>:${CVODE_INCLUDE_DIRS}>
)

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)
if(USE_LIS)
    target_compile_definitions(MathLib PUBLIC USE_LIS)
endif()
if(OGS_USE_MKL)
    target_compile_definitions(MathLib PRIVATE USE_MKL)
endif()
if(OGS_USE_PETSC)
    target_compile_definitions(MathLib PUBLIC USE_PETSC)
endif()
if(OGS_USE_EIGEN)
    target_compile_definitions(MathLib PUBLIC OGS_USE_EIGEN)
endif()

target_link_libraries(MathLib
    PUBLIC
        BaseLib
        ${OpenMP_CXX_LIBRARIES}
        $<$<TARGET_EXISTS:LAPACK::LAPACK>:LAPACK::LAPACK>
        $<$<BOOL:${OGS_USE_LIS}>:${LIS_LIBRARIES}>
        $<$<BOOL:${OGS_USE_MKL}>:${MKL_LIBRARIES}>
    PRIVATE
        spdlog::spdlog
)

if(OGS_USE_PETSC)
    target_link_libraries(MathLib PUBLIC ${PETSC_LIBRARIES})
endif()
if(OGS_USE_CVODE)
    target_link_libraries(MathLib PUBLIC ${CVODE_LIBRARIES})
endif()
