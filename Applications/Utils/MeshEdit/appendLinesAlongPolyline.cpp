/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshGeoToolsLib/AppendLinesAlongPolyline.h"

#include <tclap/CmdLine.h>

#include "Applications/FileIO/readGeometryFromFile.h"
#include "BaseLib/FileTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/PolylineVec.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Append line elements into a mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the mesh will be written to", true, "",
        "file name of output mesh");
    cmd.add(mesh_out);
    TCLAP::ValueArg<std::string> geoFileArg(
        "g", "geo-file",
        "the name of the geometry file which contains polylines", true, "",
        "the name of the geometry file");
    cmd.add(geoFileArg);

    TCLAP::ValueArg<std::string> gmsh_path_arg("g", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(gmsh_path_arg);

    // parse arguments
    cmd.parse(argc, argv);

    // read GEO objects
    GeoLib::GEOObjects geo_objs;
    FileIO::readGeometryFromFile(geoFileArg.getValue(), geo_objs,
                                 gmsh_path_arg.getValue());

    auto const geo_names = geo_objs.getGeometryNames();
    if (geo_names.empty())
    {
        ERR("No geometries found.");
        return EXIT_FAILURE;
    }
    const GeoLib::PolylineVec* ply_vec(
        geo_objs.getPolylineVecObj(geo_names[0]));
    if (!ply_vec)
    {
        ERR("Could not find polylines in geometry '{:s}'.", geo_names.front());
        return EXIT_FAILURE;
    }

    // read a mesh
    MeshLib::Mesh const* const mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    if (!mesh)
    {
        ERR("Mesh file '{:s}' not found", mesh_in.getValue());
        return EXIT_FAILURE;
    }
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    // add line elements
    std::unique_ptr<MeshLib::Mesh> new_mesh =
        MeshGeoToolsLib::appendLinesAlongPolylines(*mesh, *ply_vec);
    INFO("Mesh created: {:d} nodes, {:d} elements.",
         new_mesh->getNumberOfNodes(), new_mesh->getNumberOfElements());

    MeshLib::IO::writeMeshToFile(*new_mesh, mesh_out.getValue());

    return EXIT_SUCCESS;
}
