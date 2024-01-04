/**
 * \file
 * \brief Extracts the entire boundary from the given mesh.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshSurfaceExtraction.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Tool extracts the boundary of the given mesh. The documentation is "
        "available at "
        "https://docs.opengeosys.org/docs/tools/meshing-submeshes/"
        "extract-boundary.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the surface mesh should be written to", false, "",
        "file name of output mesh");
    cmd.add(mesh_out);

    TCLAP::SwitchArg use_ascii_arg("", "ascii-output",
                                   "If the switch is set use ascii instead of "
                                   "binary format for data in the vtu output.",
                                   false);
    cmd.add(use_ascii_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    std::unique_ptr<MeshLib::Mesh const> mesh(MeshLib::IO::readMeshFromFile(
        mesh_in.getValue(), true /* compute_element_neighbors */));

    if (!mesh)
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    // extract surface
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshToolsLib::BoundaryExtraction::getBoundaryElementsAsMesh(
            *mesh,
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Face)));

    INFO("Created surface mesh: {:d} nodes, {:d} elements.",
         surface_mesh->getNumberOfNodes(), surface_mesh->getNumberOfElements());

    std::string out_fname(mesh_out.getValue());
    if (out_fname.empty())
    {
        out_fname = BaseLib::dropFileExtension(mesh_in.getValue()) + "_sfc.vtu";
    }

    auto const data_mode =
        use_ascii_arg.getValue() ? vtkXMLWriter::Ascii : vtkXMLWriter::Binary;
    MeshLib::IO::writeVtu(*surface_mesh, out_fname, data_mode);

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
