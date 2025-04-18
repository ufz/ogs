# Source files grouped by a directory
get_source_files(SOURCES)
append_source_files(SOURCES PhreeqcIOData)
append_source_files(SOURCES PhreeqcKernelData)
append_source_files(SOURCES SelfContainedSolverData)
append_source_files(SOURCES Common)

# Create the library
ogs_add_library(ChemistryLib ${SOURCES})

target_link_libraries(
    ChemistryLib PUBLIC Eigen3::Eigen IPhreeqc PRIVATE BaseLib MaterialLib
                                                       MathLib MeshLib NumLib
)

# See https://github.com/ufz/ogs/pull/2982#issuecomment-641086788
set_source_files_properties(
    PhreeqcIO.cpp CreateChemicalSolverInterface.cpp PhreeqcKernel.cpp
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE
)

target_precompile_headers(
    ChemistryLib PRIVATE [["BaseLib/Error.h"]] [["BaseLib/ConfigTree.h"]]
    [["BaseLib/Logging.h"]]
)
