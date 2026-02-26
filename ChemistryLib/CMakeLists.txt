# Source files grouped by a directory
get_source_files(SOURCES)
append_source_files(SOURCES PhreeqcIOData)
append_source_files(SOURCES PhreeqcKernelData)
append_source_files(SOURCES SelfContainedSolverData)
append_source_files(SOURCES Common)

# Create the library
ogs_add_library(ChemistryLib ${SOURCES})

target_link_libraries(
    ChemistryLib PUBLIC Eigen3::Eigen IPhreeqc MaterialLib
                 PRIVATE BaseLib MathLib MeshLib NumLib
)

target_precompile_headers(
    ChemistryLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
)
