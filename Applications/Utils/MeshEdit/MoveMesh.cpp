/**
 * @file MoveMesh.cpp
 * @date Jan 17, 2014
 * @brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "GeoLib/AABB.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEditing/moveMeshNodes.h"
#include "MeshLib/Mesh.h"

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Moves the mesh nodes using the given displacement vector or if no displacement vector is given, moves the mesh nodes such that the centroid of the given mesh is in the origin.", ' ', "0.1");
    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-m meshfile".
    TCLAP::ValueArg<std::string> mesh_arg("m","mesh","input mesh file",true,"","string");

    // Add the argument mesh_arg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add( mesh_arg );

    TCLAP::ValueArg<double> x_arg("x","x","displacement in x direction", false, 0.0,"floating point number");
    cmd.add(x_arg);
    TCLAP::ValueArg<double> y_arg("y","y","displacement in y direction", false, 0.0,"floating point number");
    cmd.add(y_arg);
    TCLAP::ValueArg<double> z_arg("z","z","displacement in z direction", false, 0.0,"floating point number");
    cmd.add(z_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg("o","output-mesh","output mesh file", false, "", "string");
    cmd.add(mesh_out_arg);

    cmd.parse( argc, argv );

    std::string fname (mesh_arg.getValue());

    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(fname));

    if (!mesh) {
        ERR("Could not read mesh from file \"%s\".", fname.c_str());
        return EXIT_FAILURE;
    }

    MeshLib::Node displacement(0.0, 0.0, 0.0);
    if (fabs(x_arg.getValue()) < std::numeric_limits<double>::epsilon()
        && fabs(y_arg.getValue()) < std::numeric_limits<double>::epsilon()
        && fabs(z_arg.getValue()) < std::numeric_limits<double>::epsilon()) {
        GeoLib::AABB aabb(mesh->getNodes().begin(), mesh->getNodes().end());
        displacement[0] = -(aabb.getMaxPoint()[0] + aabb.getMinPoint()[0])/2.0;
        displacement[1] = -(aabb.getMaxPoint()[1] + aabb.getMinPoint()[1])/2.0;
        displacement[2] = -(aabb.getMaxPoint()[2] + aabb.getMinPoint()[2])/2.0;
    } else {
        displacement[0] = x_arg.getValue();
        displacement[1] = y_arg.getValue();
        displacement[2] = z_arg.getValue();
    }

    INFO("translate model (%f, %f, %f).",
         displacement[0],
         displacement[1],
         displacement[2]);
    MeshLib::moveMeshNodes(
        mesh->getNodes().begin(), mesh->getNodes().end(), displacement);

    std::string out_fname(mesh_out_arg.getValue());
    if (out_fname.empty()) {
        out_fname = BaseLib::dropFileExtension(mesh_out_arg.getValue());
        out_fname += "_displaced.vtu";
    }

    MeshLib::IO::writeMeshToFile(*mesh, out_fname);

    return EXIT_SUCCESS;
}
