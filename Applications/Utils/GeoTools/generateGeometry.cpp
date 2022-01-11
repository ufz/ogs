/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ThirdParty
#include <tclap/CmdLine.h>

#include <numeric>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Utils.h"
#include "InfoLib/GitInfo.h"

std::tuple<std::vector<GeoLib::Polyline*>, GeoLib::PolylineVec::NameIdMap>
appendNamedPolyline(std::unique_ptr<GeoLib::Polyline> polyline,
                    std::string&& polyline_name)
{
    std::vector<GeoLib::Polyline*> lines;
    GeoLib::PolylineVec::NameIdMap name_map;

    lines.push_back(polyline.release());
    name_map[std::move(polyline_name)] = lines.size() - 1;

    return {lines, name_map};
}

void generateSinglePointGeometry(GeoLib::Point const& point,
                                 std::string&& point_name,
                                 std::string& geometry_name,
                                 GeoLib::GEOObjects& geometry)
{
    std::vector<GeoLib::Point*> points;
    points.push_back(new GeoLib::Point{point});

    GeoLib::PointVec::NameIdMap name_map{{std::move(point_name), 0}};

    geometry.addPointVec(std::move(points), geometry_name, std::move(name_map));
}

void generatePolylineGeometry(GeoLib::Point const& point0,
                              GeoLib::Point const& point1,
                              int const number_of_subdivisions,
                              std::string&& polyline_name,
                              std::string& geometry_name,
                              GeoLib::GEOObjects& geometry)
{
    auto intermediate_points = GeoLib::generateEquidistantPoints(
        point0, point1, number_of_subdivisions);
    std::vector<GeoLib::Point*> points(intermediate_points.begin(),
                                       intermediate_points.end());
    geometry.addPointVec(std::move(points), geometry_name,
                         GeoLib::PointVec::NameIdMap{});
    auto const& point_vec = *geometry.getPointVecObj(geometry_name);

    std::vector<std::size_t> polyline_point_ids(point_vec.getVector().size());
    std::iota(begin(polyline_point_ids), end(polyline_point_ids), 0);

    auto [lines, name_map] = appendNamedPolyline(
        GeoLib::createPolyline(point_vec, std::move(polyline_point_ids)),
        std::move(polyline_name));

    geometry.addPolylineVec(std::move(lines), geometry_name,
                            std::move(name_map));
}

std::vector<GeoLib::Point*> generateQuadPoints(
    std::array<GeoLib::Point, 4> const& points,
    std::array<int, 4> const& number_of_subdivisions_per_edge)
{
    std::vector<GeoLib::Point*> quad_points;

    auto addPointsOnLine = [&quad_points](auto const& begin, auto const& end,
                                          auto const number_of_subdivisions)
    {
        auto intermediate_points = GeoLib::generateEquidistantPoints(
            begin, end, number_of_subdivisions);
        quad_points.insert(quad_points.end(), intermediate_points.begin(),
                           --intermediate_points.end());
        delete intermediate_points.back();  // Release last point, other points
                                            // are managed by GEOObjects.
    };

    addPointsOnLine(points[0], points[1], number_of_subdivisions_per_edge[0]);
    addPointsOnLine(points[1], points[2], number_of_subdivisions_per_edge[1]);
    addPointsOnLine(points[2], points[3], number_of_subdivisions_per_edge[2]);
    addPointsOnLine(points[3], points[0], number_of_subdivisions_per_edge[3]);

    return quad_points;
}

int generateQuadGeometry(GeoLib::Point const& point0,
                         GeoLib::Point const& point1,
                         int const number_of_subdivisions_first_x,
                         int const number_of_subdivisions_second_x,
                         int const number_of_subdivisions_first_y,
                         int const number_of_subdivisions_second_y,
                         int const number_of_subdivisions_first_z,
                         int const number_of_subdivisions_second_z,
                         std::string&& quad_name, std::string& geometry_name,
                         GeoLib::GEOObjects& geometry)
{
    std::array<GeoLib::Point, 4> edge_points;
    edge_points[0] = point0;
    edge_points[2] = point1;
    std::array<int, 4> number_of_subdivisions;
    if (point0[0] != point1[0] && point0[1] != point1[1] &&
        point0[2] == point1[2])
    {
        // quad in xy plane
        edge_points[1] =
            GeoLib::Point{point1[0], point0[1], point0[2]};  // right front
        edge_points[3] =
            GeoLib::Point{point0[0], point1[1], point0[2]};  // left back
        number_of_subdivisions = {
            number_of_subdivisions_first_x, number_of_subdivisions_first_y,
            number_of_subdivisions_second_x, number_of_subdivisions_second_y};
    }
    else if (point0[0] != point1[0] && point0[1] == point1[1] &&
             point0[2] != point1[2])
    {
        // quad in xz plane
        edge_points[1] =
            GeoLib::Point{point1[0], point1[1], point0[2]};  // lower right
        edge_points[3] =
            GeoLib::Point{point0[0], point0[1], point1[2]};  // upper left
        number_of_subdivisions = {
            number_of_subdivisions_first_x, number_of_subdivisions_first_z,
            number_of_subdivisions_second_x, number_of_subdivisions_second_z};
    }
    else if (point0[0] == point1[0] && point0[1] != point1[1] &&
             point0[2] != point1[2])
    {
        // quad in yz plane
        edge_points[1] =
            GeoLib::Point{point1[0], point1[1], point0[2]};  // lower back
        edge_points[3] =
            GeoLib::Point{point0[0], point0[1], point1[2]};  // upper front
        number_of_subdivisions = {
            number_of_subdivisions_first_y, number_of_subdivisions_first_z,
            number_of_subdivisions_second_y, number_of_subdivisions_second_z};
    }
    else
    {
        ERR("Input coordinates don't describe an axis aligned polyline or "
            "quad.");
        return EXIT_FAILURE;
    }

    geometry.addPointVec(
        generateQuadPoints(edge_points, number_of_subdivisions), geometry_name,
        GeoLib::PointVec::NameIdMap{});
    auto const& point_vec = *geometry.getPointVecObj(geometry_name);

    std::vector<std::size_t> polyline_point_ids(point_vec.getVector().size() +
                                                1);
    std::iota(begin(polyline_point_ids), end(polyline_point_ids), 0);
    polyline_point_ids.back() = polyline_point_ids.front();

    auto [lines, name_map] = appendNamedPolyline(
        GeoLib::createPolyline(point_vec, std::move(polyline_point_ids)),
        std::move(quad_name));

    geometry.addPolylineVec(std::move(lines), geometry_name,
                            std::move(name_map));
    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Generate point, axis parallel polyline or axis parallel quad "
        "geometry. \n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<double> z0("", "z0", "z coordinate of the first point",
                               true, 0.0, "z0");
    cmd.add(z0);
    TCLAP::ValueArg<double> y0("", "y0", "y coordinate of the first point",
                               true, 0.0, "y0");
    cmd.add(y0);
    TCLAP::ValueArg<double> x0("", "x0", "x coordinate of the first point",
                               true, 0.0, "x0");
    cmd.add(x0);
    TCLAP::ValueArg<double> z1("", "z1", "z coordinate of the first point",
                               true, 1.0, "z1");
    cmd.add(z1);
    TCLAP::ValueArg<double> y1("", "y1", "y coordinate of the first point",
                               true, 1.0, "y1");
    cmd.add(y1);
    TCLAP::ValueArg<double> x1("", "x1", "x coordinate of the first point",
                               true, 1.0, "x1");
    cmd.add(x1);
    TCLAP::ValueArg<int> nx("", "nx", "number of subdivisions in x direction",
                            false, 0, "number of subdivisions in x direction");
    cmd.add(nx);
    TCLAP::ValueArg<int> nx1("", "nx1", "number of subdivisions in x direction",
                             false, -1,
                             "number of subdivisions in x direction");
    cmd.add(nx1);
    TCLAP::ValueArg<int> ny("", "ny", "number of subdivisions in y direction",
                            false, 0, "number of subdivisions in y direction");
    cmd.add(ny);
    TCLAP::ValueArg<int> ny1("", "ny1", "number of subdivisions in y direction",
                             false, -1,
                             "number of subdivisions in y direction");
    cmd.add(ny1);
    TCLAP::ValueArg<int> nz("", "nz", "number of subdivisions in z direction",
                            false, 0, "number of subdivisions in z direction");
    cmd.add(nz);
    TCLAP::ValueArg<int> nz1("", "nz1", "number of subdivisions in z direction",
                             false, -1,
                             "number of subdivisions in z direction");
    cmd.add(nz1);
    TCLAP::ValueArg<std::string> geometry_name(
        "", "geometry_name", "name of the generated geometry", false,
        "conceptual model", "name of the geometry");
    cmd.add(geometry_name);
    TCLAP::ValueArg<std::string> polyline_name(
        "", "polyline_name", "name of the generated polyline", false,
        "polyline", "name of the generated polyline");
    cmd.add(polyline_name);
    TCLAP::ValueArg<std::string> geo_output_arg(
        "o", "output", "output geometry file (*.gml)", true, "", "output file");
    cmd.add(geo_output_arg);
    cmd.parse(argc, argv);

    auto const p0 = GeoLib::Point{x0.getValue(), y0.getValue(), z0.getValue()};
    auto const p1 = GeoLib::Point{x1.getValue(), y1.getValue(), z1.getValue()};

    GeoLib::GEOObjects geometry;
    auto constexpr eps = std::numeric_limits<double>::epsilon();
    if (p1[0] - p0[0] < eps && p1[1] - p0[1] < eps && p1[2] - p0[2] < eps)
    {
        generateSinglePointGeometry(p0, std::move(polyline_name.getValue()),
                                    geometry_name.getValue(), geometry);
    }
    else if ((p1[0] - p0[0] >= eps && p1[1] - p0[1] < eps &&
              p1[2] - p0[2] < eps) ||
             (p1[0] - p0[0] < eps && p1[1] - p0[1] >= eps &&
              p1[2] - p0[2] < eps) ||
             (p1[0] - p0[0] < eps && p1[1] - p0[1] < eps &&
              p1[2] - p0[2] >= eps))
    {
        generatePolylineGeometry(p0, p1, nx.getValue(),
                                 std::move(polyline_name.getValue()),
                                 geometry_name.getValue(), geometry);
    }
    else
    {
        auto eval = [](int v, int v1)
        {
            if (v1 == -1)
            {
                return v;
            }
            else
            {
                return v1;
            }
        };
        if (generateQuadGeometry(
                p0, p1, nx.getValue(), eval(nx.getValue(), nx1.getValue()),
                ny.getValue(), eval(ny.getValue(), ny1.getValue()),
                nz.getValue(), eval(nz.getValue(), nz1.getValue()),
                std::move(polyline_name.getValue()), geometry_name.getValue(),
                geometry) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }
    }

    GeoLib::IO::BoostXmlGmlInterface xml{geometry};
    xml.export_name = geometry.getGeometryNames()[0];
    BaseLib::IO::writeStringToFile(xml.writeToString(),
                                   geo_output_arg.getValue());

    return EXIT_SUCCESS;
}
