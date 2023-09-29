# Source files
get_source_files(SOURCES)

append_source_files(SOURCES MeshEditing)
append_source_files(SOURCES MeshGenerators)
append_source_files(SOURCES MeshQuality)

# Create the library
ogs_add_library(MeshToolsLib ${SOURCES})

target_link_libraries(
    MeshToolsLib PUBLIC GeoLib
                        MeshLib
                        range-v3
                        VTK::IOXML
                        VTK::IOLegacy
                        $<$<TARGET_EXISTS:VTK::IOParallelXML>:VTK::IOParallelXML>
                        $<$<TARGET_EXISTS:VTK::ParallelMPI>:VTK::ParallelMPI>
                 PRIVATE BaseLib MathLib NumLib GitInfoLib
)

target_precompile_headers(MeshToolsLib PRIVATE [["BaseLib/Error.h"]]
    [["BaseLib/ConfigTree.h"]] [["BaseLib/Logging.h"]] [["MeshLib/Mesh.h"]]
    [["MeshLib/Elements/Element.h"]])
