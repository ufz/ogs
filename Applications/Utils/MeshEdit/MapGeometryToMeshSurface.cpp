
/*
 * \date 2015-04-20
 * \brief Map geometric objects to the surface of the given mesh.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 */

#include <algorithm>
#include <cstdlib>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/IO/readMeshFromFile.h"

#include "MeshGeoToolsLib/GeoMapper.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Maps geometric objects to the surface of a given mesh."
        "The documentation is available at https://docs.opengeosys.org/docs/tools/model-preparation/map-geometric-object-to-the-surface-of-a-mesh",
        ' ',
        "0.1");
    TCLAP::ValueArg<std::string> mesh_in("m", "mesh-file",
        "the name of the file containing the mesh", true,
        "", "file name");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> input_geometry_fname("i", "input-geometry",
        "the name of the file containing the input geometry", true,
        "", "file name");
    cmd.add(input_geometry_fname);
    TCLAP::ValueArg<bool> additional_insert_mapping("a", "additional-insert-mapping",
        "if true advanced mapping algorithm will be applied, i.e. a new "
        "geometry will be created and possibly new points will be inserted.", false,
        true, "boolean value");
    cmd.add(additional_insert_mapping);
    TCLAP::ValueArg<std::string> output_geometry_fname("o", "output-geometry",
        "the name of the file containing the input geometry", true,
        "", "file name");
    cmd.add(output_geometry_fname);
    cmd.parse(argc, argv);

    // *** read geometry
    GeoLib::GEOObjects geometries;
    {
        GeoLib::IO::BoostXmlGmlInterface xml_io(geometries);
        if (xml_io.readFile(input_geometry_fname.getValue())) {
            INFO("Read geometry from file \"%s\".",
                input_geometry_fname.getValue().c_str());
        } else {
            return EXIT_FAILURE;
        }
    }

    std::string geo_name;
    {
        std::vector<std::string> geo_names;
        geometries.getGeometryNames(geo_names);
        geo_name = geo_names[0];
    }

    MeshGeoToolsLib::GeoMapper geo_mapper(geometries, geo_name);

    // *** read mesh
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));

    if (additional_insert_mapping.getValue()) {
        geo_mapper.advancedMapOnMesh(*mesh);
    } else {
        geo_mapper.mapOnMesh(mesh.get());
    }

    {
        GeoLib::IO::BoostXmlGmlInterface xml_io(geometries);
        xml_io.setNameForExport(geo_name);
        xml_io.writeToFile(output_geometry_fname.getValue());
    }
    return EXIT_SUCCESS;
}
