/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include <memory>
#include <gtest/gtest.h>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/GeoMapper.h"
#include "Tests/AutoCheckTools.h"

namespace ac = autocheck;

struct MeshGeoToolsLibGeoMapper : public ::testing::Test
{
    ac::IntervalGenerator<double> x_gen{0, 1};
    ac::IntervalGenerator<double> y_gen{0, 1};
    ac::IntervalGenerator<double> z_gen{-10, 10};
    ac::IntervalTupleGenerator<double> tuple_gen{x_gen, y_gen, z_gen};
    ac::cons_generator<GeoLib::Point, ac::IntervalTupleGenerator<double>>
        points_gen{tuple_gen};

    // Generates structured surface mesh, approximation of the surface described
    // by the given function, i.e., std::cos(x+y).
    std::unique_ptr<MeshLib::Mesh> _surface_mesh{
        MeshLib::MeshGenerator::createSurfaceMesh(
            "Test", MathLib::Point3d{ {{0.0, 0.0, 0.0}} },
            MathLib::Point3d{ {{1.0, 1.0, 0.0}} }, {{110,60}},
                [](double x, double y) { return std::cos(x+y); })};

    ac::gtest_reporter gtest_reporter;
};

// The test maps points with random z coordinate on the surface mesh create by
// MeshLib::MeshGenerator::createSurfaceMesh that approximates a given surface
// function. If the distance of the z coordinate of the mapped point p and the
// value of the surface function f(p_x, p_y) is smaller than a tolerance it is
// assumed that the mapping is correct.
TEST_F(MeshGeoToolsLibGeoMapper, PointsOnSurfaceMesh)
{
    auto testMapPointsOnMeshSurface = [this](
        std::vector<GeoLib::Point>& pnts) -> bool
    {
        GeoLib::GEOObjects geo_obj;
        std::string geo_name("TestGeoMapperPoints");
        auto points = std::make_unique<std::vector<GeoLib::Point*>>();
        for (auto & p : pnts) {
            points->push_back(new GeoLib::Point(p));
        }
        geo_obj.addPointVec(std::move(points), geo_name);
        MeshGeoToolsLib::GeoMapper geo_mapper(geo_obj, geo_name);

        geo_mapper.mapOnMesh(_surface_mesh.get());

        auto const mapped_points(geo_obj.getPointVec(geo_name));
        double const eps(0.01);
        for (auto pnt : *mapped_points)
        {
            GeoLib::Point const& p(*pnt);
            if (0.0 <= p[0] && p[0] <= 1.0 && 0.0 <= p[1] && p[1] <= 1.0)
            {
                if (std::abs(std::cos(p[0]+p[1]) - p[2]) >= eps)
                {
                    INFO("std::cos(%f + %f) = %f, %f",
                        p[0], p[1], cos(p[0]+p[1]), p[2]);
                    return false;
                }
            }
        }
        return true;
    };

    ac::check<std::vector<GeoLib::Point>>(
        testMapPointsOnMeshSurface,
        10,
        ac::make_arbitrary(ac::fix(100,list_of(points_gen))),
        gtest_reporter);
}

