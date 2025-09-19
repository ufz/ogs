/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <array>
#include <memory>
#include <string>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/MeshRevision.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Revises meshes by collapsing duplicate nodes, (optionally) removing "
        "elements of lower dimension than the mesh itself, and subdividing "
        "geometrically broken mesh elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<unsigned> minDim_arg(
        "d", "min-ele-dim",
        "Minimum dimension of elements to be inserted into new mesh", false, 1,
        "MIN_DIMENSION");
    cmd.add(minDim_arg);

    TCLAP::ValueArg<double> eps_arg(
        "e", "eps", "Minimum distance for nodes not to be collapsed", false,
        std::numeric_limits<double>::epsilon(), "MIN_DISTANCE");
    cmd.add(eps_arg);

    TCLAP::ValueArg<std::string> output_arg("o", "output-mesh-file",
                                            "Output (.vtu) mesh file", true, "",
                                            "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg("i", "input-mesh-file",
                                           "Input (.vtu) mesh file", true, "",
                                           "INPUT_FILE");
    cmd.add(input_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    // read a mesh file
    std::unique_ptr<MeshLib::Mesh> org_mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (!org_mesh)
    {
        return EXIT_FAILURE;
    }
    INFO("Mesh read: {:d} nodes, {:d} elements.", org_mesh->getNumberOfNodes(),
         org_mesh->getNumberOfElements());

    // revise the mesh
    INFO("Simplifying the mesh...");
    MeshToolsLib::MeshRevision const rev(const_cast<MeshLib::Mesh&>(*org_mesh));
    unsigned int minDim =
        (minDim_arg.isSet() ? minDim_arg.getValue() : org_mesh->getDimension());
    std::unique_ptr<MeshLib::Mesh> new_mesh(
        rev.simplifyMesh("revised_mesh", eps_arg.getValue(), minDim));

    // write into a file
    if (new_mesh)
    {
        INFO("Revised mesh: {:d} nodes, {:d} elements.",
             new_mesh->getNumberOfNodes(), new_mesh->getNumberOfElements());
        MeshLib::IO::writeMeshToFile(*new_mesh, output_arg.getValue());
    }

    return EXIT_SUCCESS;
}
