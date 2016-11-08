/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <memory>
#include <string>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/QuadraticeMeshGenerator.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"


int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Create quadratic order mesh", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input_arg("i", "input-mesh-file","input mesh file",true,"","string");
    cmd.add( input_arg );
    TCLAP::ValueArg<std::string> output_arg("o", "output-mesh-file","output mesh file",true,"","string");
    cmd.add( output_arg );
    cmd.parse( argc, argv );

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (!mesh)
        return EXIT_FAILURE;

    INFO("Create a quadratic order mesh");
    auto new_mesh(MeshLib::createQuadraticOrderMesh(*mesh));

    INFO("Save the new mesh into a file");
    MeshLib::IO::writeMeshToFile(*new_mesh, output_arg.getValue());

    return EXIT_SUCCESS;
}
