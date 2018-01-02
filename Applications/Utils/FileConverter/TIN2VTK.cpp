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
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"

// GeoLib
#include "GeoLib/Point.h"
#include "GeoLib/Surface.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/IO/TINInterface.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/convertMeshToGeo.h"


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts TIN file into VTU file.", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> inArg("i", "input-tin-file",
                                         "the name of the file containing the input TIN", true,
                                         "", "string");
    cmd.add(inArg);
    TCLAP::ValueArg<std::string> outArg("o", "output-vtu-file",
                                          "the name of the file the mesh will be written to", true,
                                          "", "string");
    cmd.add(outArg);
    cmd.parse(argc, argv);

    INFO("reading the TIN file...");
    const std::string tinFileName(inArg.getValue());
    auto pnt_vec = std::make_unique<std::vector<GeoLib::Point*>>();
    GeoLib::PointVec point_vec("SurfacePoints", std::move(pnt_vec));
    std::unique_ptr<GeoLib::Surface> sfc(
        GeoLib::IO::TINInterface::readTIN(tinFileName, point_vec));
    if (!sfc)
        return EXIT_FAILURE;
    INFO("TIN read:  %d points, %d triangles", pnt_vec->size(), sfc->getNumberOfTriangles());

    INFO("converting to mesh data");
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::convertSurfaceToMesh(*sfc, BaseLib::extractBaseNameWithoutExtension(tinFileName), std::numeric_limits<double>::epsilon()));
    INFO("Mesh created: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    INFO("Write it into VTU");
    MeshLib::IO::VtuInterface writer(mesh.get());
    writer.writeToFile(outArg.getValue());

    return EXIT_SUCCESS;
}
