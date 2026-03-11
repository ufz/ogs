# Source files grouped by a directory
append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    Curve
    InterpolationAlgorithms
    Integration
    LinAlg
    Nonlinear
    LinAlg/Eigen
)

if(OGS_USE_LIS)
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES LinAlg/Lis LinAlg/EigenLis
    )
endif()

if(OGS_USE_PETSC)
    append_source_and_header_files(
        _sources _headers _private DIRECTORIES LinAlg/PETSc
    )
endif()

# Create the library
ogs_add_library(
    MathLib
    GENERATE_EXPORT_HEADER
    ${_sources}
    PUBLIC_HEADERS
    ${_headers}
    PRIVATE_HEADERS
    KelvinVector-impl.h
)

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(
    MathLib
    PUBLIC
        BaseLib
        $<$<TARGET_EXISTS:LAPACK::LAPACK>:LAPACK::LAPACK>
        $<$<BOOL:${OGS_USE_LIS}>:${LIS_LIBRARIES}>
        $<$<BOOL:${OGS_USE_PETSC}>:PkgConfig::PETSC>
        Eigen3::Eigen
        Boost::math
        $<$<AND:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,$<BOOL:$<STREQUAL:${OGS_EIGEN_PARALLEL_BACKEND},OpenMP>>>:OpenMP::OpenMP_CXX>
    PRIVATE $<$<BOOL:${OGS_USE_MKL}>:MKL::MKL>
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
        $<$<BOOL:${EIGEN_NO_DEBUG}>:EIGEN_NO_DEBUG>
        $<$<BOOL:${EIGEN_DONT_VECTORIZE}>:EIGEN_DONT_VECTORIZE>
        $<$<IF:$<STREQUAL:${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES},Default>,$<NOT:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>>,$<BOOL:${OGS_EIGEN_DYNAMIC_SHAPE_MATRICES}>>:OGS_EIGEN_DYNAMIC_SHAPE_MATRICES>
        $<$<BOOL:${OGS_USE_EIGEN_UNSUPPORTED}>:USE_EIGEN_UNSUPPORTED>
        $<$<BOOL:${OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN}>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
        $<$<AND:$<BOOL:${OGS_USE_MKL}>,$<BOOL:$<STREQUAL:${OGS_EIGEN_PARALLEL_BACKEND},MKL>>>:EIGEN_USE_MKL_ALL>
)
if(OGS_USE_MKL)
    set_source_files_properties(
        LinAlg/Eigen/EigenLinearSolver.cpp PROPERTIES COMPILE_DEFINITIONS
                                                      USE_MKL
    )
endif()

target_precompile_headers(
    MathLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]] <Eigen/Core>
)
