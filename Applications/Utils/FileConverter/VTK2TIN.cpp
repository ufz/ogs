/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <memory>
#include <string>
#include <fstream>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"

// GeoLib
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Surface.h"
#include "GeoLib/IO/TINInterface.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/convertMeshToGeo.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts VTK mesh into TIN file.", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
                                         "the name of the file containing the input mesh", true,
                                         "", "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out("o", "TIN-output-file",
                                          "the name of the file the TIN will be written to", true,
                                          "", "file name of output TIN");
    cmd.add(mesh_out);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::IO::VtuInterface::readVTUFile(mesh_in.getValue()));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    INFO("Converting the mesh to TIN");
    GeoLib::GEOObjects geo_objects;
    if (MeshLib::convertMeshToGeo(*mesh, geo_objects)) {
        INFO("Writing TIN into the file");
        GeoLib::IO::TINInterface::writeSurfaceAsTIN(
            *(*geo_objects.getSurfaceVec(mesh->getName()))[0],
            mesh_out.getValue());
    }

    return EXIT_SUCCESS;
}
