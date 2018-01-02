/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <array>
#include <memory>
#include <string>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/MeshRevision.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Mesh revision tool", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input_arg("i", "input-mesh-file","input mesh file",true,"","string");
    cmd.add( input_arg );
    TCLAP::ValueArg<std::string> output_arg("o", "output-mesh-file","output mesh file",true,"","string");
    cmd.add( output_arg );
    TCLAP::SwitchArg simplify_arg("s","simplify","simplify the mesh (removing duplicated nodes)");
    cmd.add( simplify_arg );
    TCLAP::ValueArg<double> eps_arg("e", "eps","Minimum distance for nodes not to be collapsed",
                                    false, std::numeric_limits<double>::epsilon(),"float");
    cmd.add( eps_arg );
    TCLAP::ValueArg<unsigned> minDim_arg("d", "min-ele-dim","Minimum dimension of elements to be inserted into new mesh",
                                    false, 1, "unsigned");
    cmd.add( minDim_arg );
    cmd.parse( argc, argv );

    // read a mesh file
    std::unique_ptr<MeshLib::Mesh> org_mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (!org_mesh)
        return EXIT_FAILURE;
    INFO("Mesh read: %d nodes, %d elements.", org_mesh->getNumberOfNodes(), org_mesh->getNumberOfElements());

    // revise the mesh
    std::unique_ptr<MeshLib::Mesh> new_mesh;
    if (simplify_arg.getValue()) {
        INFO("Simplifying the mesh...");
        MeshLib::MeshRevision rev(const_cast<MeshLib::Mesh&>(*org_mesh));
        unsigned int minDim = (minDim_arg.isSet() ? minDim_arg.getValue() : org_mesh->getDimension());
        new_mesh.reset(
            rev.simplifyMesh("revised_mesh", eps_arg.getValue(), minDim));
    }

    // write into a file
    if (new_mesh) {
        INFO("Revised mesh: %d nodes, %d elements.", new_mesh->getNumberOfNodes(), new_mesh->getNumberOfElements());
        MeshLib::IO::writeMeshToFile(*new_mesh, output_arg.getValue());
    }

    return EXIT_SUCCESS;
}
