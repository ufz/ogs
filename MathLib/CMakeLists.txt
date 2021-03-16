# Source files grouped by a directory
get_source_files(SOURCES)

append_source_files(SOURCES Curve)
append_source_files(SOURCES InterpolationAlgorithms)
append_source_files(SOURCES Integration)
append_source_files(SOURCES LinAlg)
if(OGS_USE_CVODE)
    append_source_files(SOURCES ODE)
endif()
append_source_files(SOURCES Nonlinear)
append_source_files(SOURCES LinAlg/Eigen)

if(OGS_USE_LIS)
    append_source_files(SOURCES LinAlg/Lis)
    append_source_files(SOURCES LinAlg/EigenLis)
endif()

if(OGS_USE_PETSC)
    append_source_files(SOURCES LinAlg/PETSc)
endif()

# Create the library
ogs_add_library(MathLib ${SOURCES})
target_include_directories(
    MathLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                   $<$<BOOL:${OGS_USE_PETSC}>:${PETSC_INCLUDE_DIRS}>
)

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(
    MathLib
    PUBLIC BaseLib
           ${OpenMP_CXX_LIBRARIES}
           $<$<TARGET_EXISTS:LAPACK::LAPACK>:LAPACK::LAPACK>
           $<$<BOOL:${OGS_USE_LIS}>:${LIS_LIBRARIES}>
           $<$<BOOL:${OGS_USE_CVODE}>:CVODE::CVODE>
           Eigen3::Eigen
    PRIVATE spdlog::spdlog
)

if(OGS_USE_PETSC)
    target_link_libraries(MathLib PUBLIC ${PETSC_LIBRARIES})
endif()
if(OGS_USE_MKL)
    target_link_libraries(MathLib PUBLIC ${MKL_LIBRARIES})
endif()

target_compile_definitions(
    MathLib
    PUBLIC $<$<BOOL:${USE_LIS}>:USE_LIS> $<$<BOOL:${OGS_USE_PETSC}>:USE_PETSC>
           $<$<BOOL:${OGS_USE_CVODE}>:CVODE_FOUND>
    PRIVATE $<$<BOOL:${OGS_USE_MKL}>:USE_MKL>
)
