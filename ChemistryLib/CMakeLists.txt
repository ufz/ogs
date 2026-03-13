# Source files grouped by a directory
append_source_and_header_files(
    _sources
    _headers
    _private
    DIRECTORIES
    .
    PhreeqcIOData
    PhreeqcKernelData
    SelfContainedSolverData
    Common
)

# Create the library
ogs_add_library(ChemistryLib ${_sources} PUBLIC_HEADERS ${_headers})

target_link_libraries(
    ChemistryLib PUBLIC Eigen3::Eigen IPhreeqc MaterialLib
    PRIVATE BaseLib MathLib MeshLib NumLib
)

target_precompile_headers(
    ChemistryLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
)
