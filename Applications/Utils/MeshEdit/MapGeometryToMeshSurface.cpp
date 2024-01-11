/*
 * \file
 * \date 2015-04-20
 * \brief Map geometric objects to the surface of the given mesh.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <algorithm>
#include <cstdlib>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Maps geometric objects to the surface of a given mesh. "
        "The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/model-preparation/"
        "map-geometric-object-to-the-surface-of-a-mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "m", "mesh-file", "the name of the file containing the mesh", true, "",
        "file name");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> input_geometry_fname(
        "i", "input-geometry",
        "the name of the file containing the input geometry", true, "",
        "file name");
    cmd.add(input_geometry_fname);
    TCLAP::SwitchArg additional_insert_mapping(
        "a", "additional-insert-mapping",
        "Advanced mapping algorithm will be applied, i.e. a new geometry will "
        "be created and possibly new points will be inserted.");
    cmd.add(additional_insert_mapping);
    TCLAP::ValueArg<std::string> output_geometry_fname(
        "o", "output-geometry",
        "the name of the file containing the input geometry", true, "",
        "file name");
    cmd.add(output_geometry_fname);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    // *** read geometry
    GeoLib::GEOObjects geometries;
    {
        GeoLib::IO::BoostXmlGmlInterface xml_io(geometries);
        if (xml_io.readFile(input_geometry_fname.getValue()))
        {
            INFO("Read geometry from file '{:s}'.",
                 input_geometry_fname.getValue());
        }
        else
        {
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return EXIT_FAILURE;
        }
    }

    auto const geo_name = geometries.getGeometryNames()[0];

    MeshGeoToolsLib::GeoMapper geo_mapper(geometries, geo_name);

    // *** read mesh
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(
        mesh_in.getValue(), true /* compute_element_neighbors */));

    if (additional_insert_mapping.getValue())
    {
        geo_mapper.advancedMapOnMesh(*mesh);
    }
    else
    {
        geo_mapper.mapOnMesh(mesh.get());
    }

    {
        GeoLib::IO::BoostXmlGmlInterface xml_io(geometries);
        xml_io.export_name = geo_name;
        BaseLib::IO::writeStringToFile(xml_io.writeToString(),
                                       output_geometry_fname.getValue());
    }
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
