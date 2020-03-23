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
target_include_directories(MathLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(MathLib PUBLIC BaseLib ${OpenMP_CXX_LIBRARIES})

if(CVODE_FOUND)
    target_include_directories(MathLib PRIVATE ${CVODE_INCLUDE_DIRS})
    target_link_libraries(MathLib PUBLIC ${CVODE_LIBRARIES})
endif()

if(LAPACK_FOUND)
    target_link_libraries(MathLib PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if(OGS_USE_LIS)
    target_link_libraries(MathLib PUBLIC ${LIS_LIBRARIES})
endif()

if(OGS_USE_MKL)
    target_link_libraries(MathLib PUBLIC ${MKL_LIBRARIES})
endif()

if(OGS_USE_PETSC)
    target_include_directories(MathLib PUBLIC ${PETSC_INCLUDE_DIRS})
    target_link_libraries(MathLib PUBLIC ${PETSC_LIBRARIES})
endif()
