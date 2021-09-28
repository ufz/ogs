/**
 * \file
 * \author Karsten Rink
 * \date   2015-02-02
 * \brief  Utility for triangulating polylines.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <string>

#include "Applications/FileIO/Legacy/createSurface.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/Polyline.h"
#include "InfoLib/GitInfo.h"

std::string output_question()
{
    WARN("Given polyline is not closed. Close polyline now?");
    WARN("Enter (Y)es for closing the line or (N)o for abort and press ENTER.");
    std::string input;
    std::cin >> input;
    return input;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Triangulates the specified polyline in the given geometry file.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "GML input file (*.gml)", true, "", "string");
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "GML output file (*.gml)", true, "", "string");
    TCLAP::ValueArg<std::string> name_arg(
        "n", "name", "Name of polyline in given file", true, "", "string");
    TCLAP::ValueArg<std::string> gmsh_path_arg("g", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(input_arg);
    cmd.add(name_arg);
    cmd.add(output_arg);
    cmd.add(gmsh_path_arg);
    cmd.parse(argc, argv);

    std::string const& file_name(input_arg.getValue());
    std::string const& polyline_name(name_arg.getValue());

    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    try
    {
        if (!xml.readFile(file_name))
        {
            ERR("Failed to load geometry file.");
            return EXIT_FAILURE;
        }
    }
    catch (std::runtime_error const& err)
    {
        ERR("Failed to read file `{:s}'.", file_name);
        ERR("{:s}", err.what());
        return EXIT_FAILURE;
    }

    auto const geo_name = geo_objects.getGeometryNames()[0];
    GeoLib::PolylineVec const* const line_vec(
        geo_objects.getPolylineVecObj(geo_name));
    GeoLib::Polyline* line = const_cast<GeoLib::Polyline*>(
        line_vec->getElementByName(polyline_name));

    // check if line exists
    if (line == nullptr)
    {
        ERR("No polyline found with name '{:s}'. Aborting...", polyline_name);
        return EXIT_FAILURE;
    }

    // check if polyline can be triangulated (i.e. closed + coplanar)
    if (!line->isCoplanar())
    {
        ERR("Polyline is not coplanar, no unambiguous triangulation possible. "
            "Aborting...");
        return EXIT_FAILURE;
    }

    if (!line->isClosed())
    {
        std::string input;
        while (input != "y" && input != "Y" && input != "n" && input != "N")
        {
            input = output_question();
        }

        if (input == "y" || input == "Y")
        {
            line->closePolyline();
            INFO("Polyline closed.");
        }
        else
        {
            return EXIT_FAILURE;
        }
    }

    INFO("Creating a surface by triangulation of the polyline ...");
    if (FileIO::createSurface(*line, geo_objects, geo_name,
                              gmsh_path_arg.getValue()))
    {
        INFO("\t done");
    }
    else
    {
        WARN("\t Creating a surface by triangulation of the polyline failed.");
    }
    GeoLib::SurfaceVec* sfc_vec(geo_objects.getSurfaceVecObj(geo_name));
    std::size_t const sfc_id = geo_objects.getSurfaceVec(geo_name)->size() - 1;
    std::string const surface_name(polyline_name + "_surface");
    for (std::size_t i = 1;; ++i)
    {
        std::string const new_surface_name =
            (i > 1) ? (surface_name + std::to_string(i)) : surface_name;
        if (sfc_vec->getElementByName(new_surface_name) == nullptr)
        {
            sfc_vec->setNameForElement(sfc_id, new_surface_name);
            break;
        }
    }

    // write new file
    xml.export_name = geo_name;
    BaseLib::IO::writeStringToFile(xml.writeToString(), output_arg.getValue());
    INFO("...done.");

    return EXIT_SUCCESS;
}
