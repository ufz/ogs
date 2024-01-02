/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <memory>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"
#include "InfoLib/GitInfo.h"

auto createMapping(std::vector<bool> const& used_points)
{
    auto getCountOrMax = [cnt = 0u](bool used) mutable
    { return used ? cnt++ : std::numeric_limits<std::size_t>::max(); };

    return used_points | ranges::views::transform(getCountOrMax) |
           ranges::to<std::vector>;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Remove points from geometry that are not used in any polyline nor "
        "surface. Attention: Names of points are not handled correctly in "
        "every case.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
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

    auto mark = [&](auto const* const object)
    { GeoLib::markUsedPoints(*object, used_points); };

    if (polylines)
    {
        ranges::for_each(*polylines, mark);
    }
    if (surfaces)
    {
        ranges::for_each(*surfaces, mark);
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

    auto const mapping = createMapping(used_points);

    auto reset = [&mapping](auto* object)
    { GeoLib::resetPointIDs(*object, mapping); };

    // reset point ids is polylines
    if (polylines)
    {
        ranges::for_each(*polylines, reset);
    }
    // reset point ids is surfaces
    if (surfaces)
    {
        ranges::for_each(*surfaces, reset);
    }

    std::stringstream unused_points_info;
    unused_points_info << std::fixed;

    // cleanup unused points
    for (GeoLib::Point*& point :
         ranges::views::zip(*points, used_points) |
             ranges::views::filter([](auto&& pair) { return !pair.second; }) |
             ranges::views::keys)
    {
        unused_points_info << point->getID() << " " << *point << '\n';
        delete point;
        point = nullptr;
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
