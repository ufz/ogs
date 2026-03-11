append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
        .
        Assembler
        DOF
        Fem
        Fem/CoordinatesMapping
        Fem/FiniteElement
        Fem/Integration
        Fem/ShapeFunction
        NumericalStability
        TimeStepping
        TimeStepping/Algorithms
        StaggeredCoupling
        ODESolver
        Extrapolation
)

# Create the library
ogs_add_library(
    NumLib
    GENERATE_EXPORT_HEADER
    ${_sources}
    PUBLIC_HEADERS
    ${_headers}
    PRIVATE_HEADERS
    ${_private}
)
set_target_properties(NumLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(
    NumLib PUBLIC BaseLib GeoLib MathLib MeshLib
                  $<$<TARGET_EXISTS:petsc>:petsc>
)

target_precompile_headers(
    NumLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]] <Eigen/Core>
)
