/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <memory>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Remove points from geometry that are not used in any polyline nor "
        "surface. Attention: Names of points are not handled correctly in "
        "every case.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2023, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> geo_output_arg(
        "o", "output", "output geometry", true, "", "gml file");
    cmd.add(geo_output_arg);
    TCLAP::ValueArg<std::string> geo_input_arg(
        "i", "input", "input geometry", true, "conceptual model", "gml file");
    cmd.add(geo_input_arg);
    cmd.parse(argc, argv);

    GeoLib::GEOObjects geometry;
    GeoLib::IO::BoostXmlGmlInterface xml{geometry};
    try
    {
        if (!xml.readFile(geo_input_arg.getValue()))
        {
            ERR("Failed to read file `{:s}'.", geo_input_arg.getValue());
            return EXIT_FAILURE;
        }
    }
    catch (std::runtime_error const& err)
    {
        ERR("Failed to read file `{:s}'.", geo_input_arg.getValue());
        ERR("{:s}", err.what());
        return EXIT_FAILURE;
    }

    auto const geo_name = geometry.getGeometryNames()[0];
    auto* points = const_cast<std::vector<GeoLib::Point*>*>(
        geometry.getPointVec(geo_name));
    auto* polylines = geometry.getPolylineVec(geo_name);
    auto* surfaces = geometry.getSurfaceVec(geo_name);

    // mark used points
    auto used_points = std::vector<bool>(points->size(), false);
    if (polylines)
    {
        for (auto const* polyline : *polylines)
        {
            for (std::size_t i = 0; i < polyline->getNumberOfPoints(); ++i)
            {
                used_points[polyline->getPointID(i)] = true;
            }
        }
    }
    if (surfaces)
    {
        for (auto const* surface : *surfaces)
        {
            for (std::size_t i = 0; i < surface->getNumberOfTriangles(); ++i)
            {
                auto const* triangle = (*surface)[i];
                used_points[(*triangle)[0]] = true;
                used_points[(*triangle)[1]] = true;
                used_points[(*triangle)[2]] = true;
            }
        }
    }

    auto const number_of_used_points = static_cast<std::size_t>(
        std::count(used_points.begin(), used_points.end(), true));
    if (number_of_used_points == 0)
    {
        INFO(
            "Geometry consists of points only. A new geometry file won't "
            "be written.");
        return EXIT_SUCCESS;
    }
    INFO(
        "{} points in this geometry file are not used in any polyline or "
        "surface",
        points->size() - number_of_used_points);

    std::vector<std::size_t> mapping(points->size(),
                                     std::numeric_limits<std::size_t>::max());
    std::size_t cnt = 0;
    for (std::size_t i = 0; i < points->size(); ++i)
    {
        if (used_points[i])
        {
            mapping[i] = cnt;
            cnt++;
        }
    }

    // reset point ids is polylines
    if (polylines)
    {
        for (auto* polyline : *polylines)
        {
            GeoLib::resetPointIDs(*polyline, mapping);
        }
    }
    // reset point ids is surfaces
    if (surfaces)
    {
        for (auto* surface : *surfaces)
        {
            for (std::size_t i = 0; i < surface->getNumberOfTriangles(); ++i)
            {
                auto* triangle = (*surface)[i];
                const_cast<std::size_t&>((*triangle)[0]) =
                    mapping[(*triangle)[0]];
                const_cast<std::size_t&>((*triangle)[1]) =
                    mapping[(*triangle)[1]];
                const_cast<std::size_t&>((*triangle)[2]) =
                    mapping[(*triangle)[2]];
            }
        }
    }

    auto const stations = createStations(*points, used_points);
    if (!stations.empty())
    {
        std::stringstream stations_string;
        stations_string << std::fixed;
        for (auto const* station : stations)
        {
            stations_string << "point_id_" << station->getName() << " "
                            << (*station)[0] << " " << (*station)[1] << " "
                            << (*station)[2] << std::endl;
        }
        INFO("stations:\n{}", stations_string.str());
    }
    BaseLib::cleanupVectorElements(stations);

    // cleanup unused points
    for (std::size_t i = 0; i < points->size(); ++i)
    {
        if (!used_points[i])
        {
            delete (*points)[i];
            (*points)[i] = nullptr;
        }
    }

    std::erase(*points, nullptr);
    if (!unused_points_info.str().empty())
    {
        INFO("Removed the following points:\n{}", unused_points_info.str());
    }

    WARN("Names of points are not handled correctly in every case.");

    xml.export_name = geometry.getGeometryNames()[0];
    BaseLib::IO::writeStringToFile(xml.writeToString(),
                                   geo_output_arg.getValue());

    return EXIT_SUCCESS;
}
