/*
 * \file
 * \brief Reset material properties in meshes in a polygonal region.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <cstdlib>
#include <vector>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "Applications/FileIO/readGeometryFromFile.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "InfoLib/GitInfo.h"
#include "MeshGeoToolsLib/MeshEditing/ResetMeshElementProperty.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshInformation.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Sets the property value of a mesh element to a given new value iff at "
        "least one node of the element is within a polygonal region that is "
        "given by a polygon. The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/model-preparation/"
        "set-properties-in-polygonal-region.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the mesh will be written to, format depends on "
        "the given file name extension",
        true, "", "file name");
    cmd.add(mesh_out);
    TCLAP::ValueArg<std::string> polygon_name_arg(
        "p", "polygon-name", "name of polygon in the geometry", true, "",
        "string");
    cmd.add(polygon_name_arg);
    TCLAP::ValueArg<std::string> geometry_fname(
        "g", "geometry",
        "the name of the file containing the input geometry (gli or gml "
        "format)",
        true, "", "file name");
    cmd.add(geometry_fname);
    TCLAP::SwitchArg any_of_arg(
        "", "any_of",
        "all nodes of an element has to be inside the polygon (default "
        "behaviour without switch) or any node of an element has to be inside "
        "(switch is given)",
        false);
    cmd.add(any_of_arg);
    TCLAP::ValueArg<char> char_property_arg(
        "c", "char-property-value", "new property value (data type char)",
        false, 'A', "character");
    cmd.add(char_property_arg);
    TCLAP::ValueArg<int> int_property_arg("i", "int-property-value",
                                          "new property value (data type int)",
                                          false, 0, "number");
    cmd.add(int_property_arg);
    TCLAP::ValueArg<bool> bool_property_arg(
        "b", "bool-property-value", "new property value (data type bool)",
        false, false, "boolean value");
    cmd.add(bool_property_arg);
    TCLAP::ValueArg<std::string> property_name_arg(
        "n", "property-name", "name of property in the mesh", false,
        "MaterialIDs", "string");
    cmd.add(property_name_arg);
    TCLAP::ValueArg<int> restrict_arg(
        "r", "restrict-to-MaterialID",
        "Restrict resetting the property to the material id", false, -1,
        "MaterialID");
    cmd.add(restrict_arg);
    TCLAP::ValueArg<std::string> mesh_in(
        "m", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> gmsh_path_arg("", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(gmsh_path_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    // *** read geometry
    GeoLib::GEOObjects geometries;
    FileIO::readGeometryFromFile(geometry_fname.getValue(), geometries,
                                 gmsh_path_arg.getValue());

    auto const geo_name = geometries.getGeometryNames()[0];

    // *** check if the data is usable
    // *** get vector of polylines
    GeoLib::PolylineVec const* plys(geometries.getPolylineVecObj(geo_name));
    if (!plys)
    {
        ERR("Could not get vector of polylines out of geometry '{:s}'.",
            geo_name);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    // *** get polygon
    GeoLib::Polyline const* ply(
        plys->getElementByName(polygon_name_arg.getValue()));
    if (!ply)
    {
        ERR("Polyline '{:s}' not found.", polygon_name_arg.getValue());
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    // *** check if the polyline is closed (i.e. is a polygon)
    if (!ply->isClosed())
    {
        ERR("Polyline '{:s}' is not closed, i.e. does not describe a region.",
            polygon_name_arg.getValue());
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    GeoLib::Polygon const polygon(*(ply));

    // *** read mesh
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    if (!mesh)
    {
        // error message written already by readMeshFromFile
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }
    std::string const& property_name(property_name_arg.getValue());

    if (char_property_arg.isSet())
    {
        MeshGeoToolsLib::resetMeshElementProperty(
            *mesh, polygon, property_name, char_property_arg.getValue(),
            restrict_arg.getValue(), any_of_arg.getValue());
    }

    if (int_property_arg.isSet())
    {
        MeshGeoToolsLib::resetMeshElementProperty(
            *mesh, polygon, property_name, int_property_arg.getValue(),
            restrict_arg.getValue(), any_of_arg.getValue());
    }

    if (bool_property_arg.isSet())
    {
        MeshGeoToolsLib::resetMeshElementProperty(
            *mesh, polygon, property_name, bool_property_arg.getValue(),
            restrict_arg.getValue(), any_of_arg.getValue());
    }

    MeshToolsLib::MeshInformation::writePropertyVectorInformation(*mesh);

    MeshLib::IO::writeMeshToFile(*mesh, mesh_out.getValue());

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
