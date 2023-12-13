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
ogs_add_library(MathLib GENERATE_EXPORT_HEADER ${SOURCES})

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(
    MathLib
    PUBLIC BaseLib
           $<$<TARGET_EXISTS:LAPACK::LAPACK>:LAPACK::LAPACK>
           $<$<BOOL:${OGS_USE_LIS}>:${LIS_LIBRARIES}>
           $<$<BOOL:${OGS_USE_CVODE}>:CVODE::CVODE>
           $<$<BOOL:${OGS_USE_PETSC}>:PkgConfig::PETSC>
           $<$<BOOL:${OGS_USE_MKL}>:MKL::MKL>
           Eigen3::Eigen
           $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
)
if(OGS_USE_MKL AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Otherwise required MKL libs get stripped out (don't show up in ldd):
    target_link_options(MathLib PUBLIC "LINKER:--no-as-needed")
endif()

if(OGS_USE_LIS)
    target_include_directories(MathLib PUBLIC ${LIS_INCLUDE_DIR})
endif()

target_compile_definitions(
    MathLib
    PUBLIC
        $<$<BOOL:${OGS_USE_LIS}>:USE_LIS>
        $<$<BOOL:${OGS_USE_CVODE}>:CVODE_FOUND>
        $<$<BOOL:${EIGEN_NO_DEBUG}>:EIGEN_NO_DEBUG>
        $<$<BOOL:${EIGEN_DONT_VECTORIZE}>:EIGEN_DONT_VECTORIZE>
        $<$<BOOL:${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_INTERNAL}>:OGS_EIGEN_DYNAMIC_SHAPE_MATRICES>
        $<$<BOOL:${OGS_USE_EIGEN_UNSUPPORTED}>:USE_EIGEN_UNSUPPORTED>
        $<$<BOOL:${OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN}>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
        $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
        $<$<AND:$<BOOL:${OGS_USE_MKL}>,$<BOOL:${OGS_EIGEN_USE_MKL}>>:EIGEN_USE_MKL_ALL>
    PRIVATE $<$<BOOL:${OGS_USE_MKL}>:USE_MKL>
)

target_precompile_headers(
    MathLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]] <Eigen/Core>
)
