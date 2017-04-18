#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES)

APPEND_SOURCE_FILES(SOURCES Curve)
APPEND_SOURCE_FILES(SOURCES InterpolationAlgorithms)
APPEND_SOURCE_FILES(SOURCES Integration)
APPEND_SOURCE_FILES(SOURCES LinAlg)
APPEND_SOURCE_FILES(SOURCES LinAlg/Dense)
APPEND_SOURCE_FILES(SOURCES LinAlg/Solvers)
APPEND_SOURCE_FILES(SOURCES ODE)
APPEND_SOURCE_FILES(SOURCES Nonlinear)

if(OGS_USE_EIGEN)
    APPEND_SOURCE_FILES(SOURCES LinAlg/Eigen)
endif()

if(OGS_USE_LIS)
    APPEND_SOURCE_FILES(SOURCES LinAlg/Lis)
    APPEND_SOURCE_FILES(SOURCES LinAlg/EigenLis)
endif()

if(OGS_USE_PETSC)
    APPEND_SOURCE_FILES(SOURCES LinAlg/PETSc)
endif()

# Create the library
add_library(MathLib ${SOURCES})

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(MathLib
    BaseLib
)

if (CVODE_FOUND)
    target_link_libraries(MathLib ${CVODE_LIBRARIES})
endif()

if(LAPACK_FOUND)
    target_link_libraries(MathLib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if (OGS_USE_LIS)
    target_link_libraries(MathLib ${LIS_LIBRARIES})
endif()

if (OGS_USE_MKL)
    target_link_libraries(MathLib ${MKL_LIBRARIES})
endif()

if (OGS_USE_PETSC)
    target_link_libraries(MathLib ${PETSC_LIBRARIES})
endif()

if(TARGET Boost)
    add_dependencies(MathLib Boost)
endif()
if(TARGET Eigen)
    add_dependencies(MathLib Eigen)
endif()

cotire(MathLib)
