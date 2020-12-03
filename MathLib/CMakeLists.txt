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
        $<$<BOOL:${CVODE_FOUND}>:${CVODE_INCLUDE_DIRS}>
)

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(MathLib
    PUBLIC
        BaseLib
        ${OpenMP_CXX_LIBRARIES}
        $<$<BOOL:${CVODE_FOUND}>:${CVODE_LIBRARIES}>
        $<$<BOOL:${LAPACK_FOUND}>:${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}>
        $<$<BOOL:${OGS_USE_LIS}>:${LIS_LIBRARIES}>
        $<$<BOOL:${OGS_USE_MKL}>:${MKL_LIBRARIES}>
        $<$<BOOL:${OGS_USE_PETSC}>:${PETSC_LIBRARIES}>
    PRIVATE
        spdlog::spdlog
)
