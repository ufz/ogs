/**
 * \file
 * \date Jan 17, 2014
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

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/moveMeshNodes.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Moves the mesh nodes using the given displacement vector or if no "
        "displacement vector is given, moves the mesh nodes such that the "
        "centroid of the given mesh is in the origin.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-m meshfile".
    TCLAP::ValueArg<std::string> mesh_arg("m", "mesh", "input mesh file", true,
                                          "", "string");

    // Add the argument mesh_arg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add(mesh_arg);

    TCLAP::ValueArg<double> x_arg("x", "x", "displacement in x direction",
                                  false, 0.0, "floating point number");
    cmd.add(x_arg);
    TCLAP::ValueArg<double> y_arg("y", "y", "displacement in y direction",
                                  false, 0.0, "floating point number");
    cmd.add(y_arg);
    TCLAP::ValueArg<double> z_arg("z", "z", "displacement in z direction",
                                  false, 0.0, "floating point number");
    cmd.add(z_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh", "output mesh file", false, "", "string");
    cmd.add(mesh_out_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif
    std::string fname(mesh_arg.getValue());

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(fname));

    if (!mesh)
    {
        ERR("Could not read mesh from file '{:s}'.", fname);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    Eigen::Vector3d displacement{x_arg.getValue(), y_arg.getValue(),
                                 z_arg.getValue()};
    if (fabs(x_arg.getValue()) < std::numeric_limits<double>::epsilon() &&
        fabs(y_arg.getValue()) < std::numeric_limits<double>::epsilon() &&
        fabs(z_arg.getValue()) < std::numeric_limits<double>::epsilon())
    {
        GeoLib::AABB aabb(mesh->getNodes().begin(), mesh->getNodes().end());
        auto const [min, max] = aabb.getMinMaxPoints();
        displacement = -(max + min) / 2.0;
    }

    INFO("translate model ({:f}, {:f}, {:f}).",
         displacement[0],
         displacement[1],
         displacement[2]);
    MeshToolsLib::moveMeshNodes(mesh->getNodes().begin(),
                                mesh->getNodes().end(), displacement);

    std::string out_fname(mesh_out_arg.getValue());
    if (out_fname.empty())
    {
        out_fname = BaseLib::dropFileExtension(mesh_out_arg.getValue());
        out_fname += "_displaced.vtu";
    }

    MeshLib::IO::writeMeshToFile(*mesh, out_fname);

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
