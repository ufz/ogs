/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

std::unique_ptr<GeoLib::GEOObjects> readGeometry(std::string const& filename)
{
    auto geo_objects = std::make_unique<GeoLib::GEOObjects>();
    GeoLib::IO::BoostXmlGmlInterface gml_reader(*geo_objects);

    DBUG("Reading geometry file '{:s}'.", filename);
    gml_reader.readFile(filename);
    return geo_objects;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts a geometry defined on a given mesh to distinct meshes. The "
        "documentation is available at "
        "https://www.opengeosys.org/docs/tools/meshing-submeshes/"
        "constructmeshesfromgeometry/.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<double> search_length_arg(
        "s",
        "searchlength",
        "search length determining radius for the node search algorithm. "
        "Non-negative floating point number (default 1e-16) ",
        false,
        1e-16,
        "float");
    cmd.add(search_length_arg);

    TCLAP::ValueArg<std::string> geometry_arg("g",
                                              "geometry",
                                              "the file name the geometry",
                                              true,
                                              "",
                                              "geometry file name");
    cmd.add(geometry_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "m",
        "mesh",
        "the file name of the mesh where the geometry is defined",
        true,
        "",
        "mesh file name");
    cmd.add(mesh_arg);

    TCLAP::SwitchArg multiple_nodes_allowed_arg(
        "", "multiple-nodes-allowed",
        "Allows multiple mesh nodes in eps environment.");
    cmd.add(multiple_nodes_allowed_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh{
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue())};

    auto const geo_objects = readGeometry(geometry_arg.getValue());

    double const search_length = search_length_arg.getValue();

    auto const extracted_meshes = constructAdditionalMeshesFromGeoObjects(
        *geo_objects,
        *mesh,
        std::make_unique<MeshGeoToolsLib::SearchLength>(search_length),
        multiple_nodes_allowed_arg.getValue());

    for (auto const& m_ptr : extracted_meshes)
    {
        if (!m_ptr)
        {
            ERR("Could not create a mesh for each given geometry.");
            return EXIT_FAILURE;
        }
        if (m_ptr->getNodes().empty())
        {
            WARN(
                "The created mesh '{:s}' hasn't any nodes or elements and thus "
                "it isn't written to file.",
                m_ptr->getName());
            continue;
        }
        MeshLib::IO::writeMeshToFile(*m_ptr, m_ptr->getName() + ".vtu");
    }

    return EXIT_SUCCESS;
}
