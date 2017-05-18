/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
#include <memory>

#include <gtest/gtest.h>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"
#include "MeshGeoToolsLib/MeshNodesAlongSurface.h"
#include "NumLib/Function/LinearInterpolationAlongPolyline.h"
#include "NumLib/Function/LinearInterpolationOnSurface.h"
#include "NumLib/Function/SpatialFunctionLinear.h"
#include "NumLib/Distribution/Distribution.h"

#include "Tests/TestTools.h"

class NumLibDistributionQuad : public testing::Test
{
public:
    NumLibDistributionQuad()
        : _geometric_size(10.0),
          _number_of_subdivisions_per_direction(10),
          _msh(MeshLib::MeshGenerator::generateRegularQuadMesh(
              _geometric_size, _number_of_subdivisions_per_direction)),
          _project_name("test"),
          _mshNodesSearcher(*_msh, MeshGeoToolsLib::SearchLength(),
                            MeshGeoToolsLib::SearchAllNodes::Yes),
          _ply0(nullptr)
    {
        // create geometry
        auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
            new std::vector<GeoLib::Point*>);
        pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
        pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
        pnts->push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));

        auto plys = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
            new std::vector<GeoLib::Polyline*>);
        _ply0 = new GeoLib::Polyline(*pnts);
        _ply0->addPoint(0);
        _ply0->addPoint(1);
        plys->push_back(_ply0);

        auto* ply1 = new GeoLib::Polyline(*pnts);
        ply1->addPoint(0);
        ply1->addPoint(1);
        ply1->addPoint(2);
        ply1->addPoint(3);
        ply1->addPoint(0);
        plys->push_back(ply1);

        auto sfcs = std::unique_ptr<std::vector<GeoLib::Surface*>>(
            new std::vector<GeoLib::Surface*>);
        _sfc1 = GeoLib::Surface::createSurface(*ply1);
        sfcs->push_back(_sfc1);

        _geo_objs.addPointVec(std::move(pnts), _project_name);
        _geo_objs.addPolylineVec(std::move(plys), _project_name);
        _geo_objs.addSurfaceVec(std::move(sfcs), _project_name);
    }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    std::unique_ptr<MeshLib::Mesh> _msh;
    GeoLib::GEOObjects _geo_objs;
    std::string _project_name;
    MeshGeoToolsLib::MeshNodeSearcher _mshNodesSearcher;
    GeoLib::Polyline* _ply0;
    GeoLib::Surface* _sfc1;
};

class NumLibDistributionHex : public testing::Test
{
public:
    NumLibDistributionHex()
        : _geometric_size(10.0),
          _number_of_subdivisions_per_direction(10),
          _msh(MeshLib::MeshGenerator::generateRegularHexMesh(
              _geometric_size, _number_of_subdivisions_per_direction)),
          _project_name("test"),
          _mshNodesSearcher(*_msh, MeshGeoToolsLib::SearchLength(),
                            MeshGeoToolsLib::SearchAllNodes::Yes),
          _ply0(nullptr)
    {
        // create geometry
        auto pnts = std::unique_ptr<std::vector<GeoLib::Point*>>(
            new std::vector<GeoLib::Point*>);
        pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
        pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
        pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
        pnts->push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
        pnts->push_back(new GeoLib::Point(0.0, 0.0, _geometric_size));
        pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, _geometric_size));
        pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, _geometric_size));
        pnts->push_back(new GeoLib::Point(0.0, _geometric_size, _geometric_size));

        auto plys = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
            new std::vector<GeoLib::Polyline*>);
        _ply0 = new GeoLib::Polyline(*pnts); // vertical polyline
        _ply0->addPoint(0);
        _ply0->addPoint(4);
        plys->push_back(_ply0);
        auto* ply1 = new GeoLib::Polyline(*pnts);  // polygon for left surface
        ply1->addPoint(0);
        ply1->addPoint(3);
        ply1->addPoint(7);
        ply1->addPoint(4);
        ply1->addPoint(0);
        plys->push_back(ply1);

        auto sfcs = std::unique_ptr<std::vector<GeoLib::Surface*>>(
            new std::vector<GeoLib::Surface*>);
        _sfc1 = GeoLib::Surface::createSurface(*ply1);
        sfcs->push_back(_sfc1);

        _geo_objs.addPointVec(std::move(pnts) ,_project_name);
        _geo_objs.addPolylineVec(std::move(plys), _project_name);
        _geo_objs.addSurfaceVec(std::move(sfcs), _project_name);
    }

protected:
    const double _geometric_size;
    const std::size_t _number_of_subdivisions_per_direction;
    std::unique_ptr<MeshLib::Mesh> _msh;
    GeoLib::GEOObjects _geo_objs;
    std::string _project_name;
    MeshGeoToolsLib::MeshNodeSearcher _mshNodesSearcher;
    GeoLib::Polyline* _ply0;
    GeoLib::Surface* _sfc1;
};

TEST_F(NumLibDistributionQuad, Linear)
{
    // f(x,y,z) = 1 + 2x + 3y + 4z
    std::array<double,4> f_coeff = {{1, 2, 3, 4}};
    MathLib::LinearFunction<double,3> linear_f(f_coeff);
    std::vector<double> expected = {{1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21}};

    NumLib::SpatialFunctionLinear f(linear_f);
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongPolyline(*_ply0).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(f, *_msh, vec_node_ids);
    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(NumLibDistributionHex, Linear)
{
    // f(x,y,z) = 1 + 2x + 3y + 4z
    std::array<double,4> f_coeff = {{1, 2, 3, 4}};
    MathLib::LinearFunction<double,3> linear_f(f_coeff);
    std::vector<double> expected(static_cast<std::size_t>(std::pow(_number_of_subdivisions_per_direction+1, 2)));
    const double dL = _geometric_size / _number_of_subdivisions_per_direction;
    for (std::size_t i=0; i<expected.size(); i++) {
        double x = 0;
        double y = (i%(_number_of_subdivisions_per_direction+1)) * dL;
        double z = (i/(_number_of_subdivisions_per_direction+1)) * dL;
        expected[i] = f_coeff[0] + f_coeff[1]*x + f_coeff[2]*y + f_coeff[3]*z;
    }

    NumLib::SpatialFunctionLinear f(linear_f);
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongSurface(*_sfc1).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(f, *_msh, vec_node_ids);

    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(NumLibDistributionQuad, InterpolationPolyline)
{
    const std::vector<std::size_t> vec_point_ids = {{0, 1}};
    const std::vector<double> vec_point_values = {{0., 100.}};
    std::vector<double> expected = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}};

    NumLib::LinearInterpolationAlongPolyline interpolate(
            *_ply0, vec_point_ids, vec_point_values, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max());
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongPolyline(*_ply0).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(interpolate, *_msh, vec_node_ids);

    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(NumLibDistributionHex, InterpolationPolyline)
{
    const std::vector<std::size_t> vec_point_ids = {{0, 4}};
    const std::vector<double> vec_point_values = {{0., 100.}};
    std::vector<double> expected = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}};

    NumLib::LinearInterpolationAlongPolyline interpolate(
            *_ply0, vec_point_ids, vec_point_values, std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::max());
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongPolyline(*_ply0).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(interpolate, *_msh, vec_node_ids);

    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(NumLibDistributionQuad, InterpolationSurface)
{
    const std::vector<std::size_t> vec_point_ids = {{0, 1, 2, 3}};
    const std::vector<double> vec_point_values = {{0., 100., 100., 0.}};
    std::vector<double> expected(_msh->getNumberOfNodes());
    for (std::size_t i=0; i<_msh->getNumberOfNodes(); i++) {
        expected[i] = static_cast<double>((i%(_number_of_subdivisions_per_direction+1)) * 10);
    }

    NumLib::LinearInterpolationOnSurface interpolate(*_sfc1, vec_point_ids, vec_point_values, std::numeric_limits<double>::max());
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongSurface(*_sfc1).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(interpolate, *_msh, vec_node_ids);

    // the machine epsilon for double is too small for this test
    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<float>::epsilon());
}

TEST_F(NumLibDistributionHex, InterpolationSurface)
{
    const std::vector<std::size_t> vec_point_ids = {{0, 3, 7, 4}};
    const std::vector<double> vec_point_values = {{0., 100., 100., 0.}};
    std::vector<double> expected(static_cast<std::size_t>(std::pow(_number_of_subdivisions_per_direction+1, 2)));
    for (std::size_t i=0; i<expected.size(); i++) {
        expected[i] = static_cast<double>((i%(_number_of_subdivisions_per_direction+1)) * 10);
    }

    NumLib::LinearInterpolationOnSurface interpolate(*_sfc1, vec_point_ids, vec_point_values, std::numeric_limits<double>::max());
    const std::vector<std::size_t>& vec_node_ids = _mshNodesSearcher.getMeshNodesAlongSurface(*_sfc1).getNodeIDs();
    std::vector<double> interpolated_values = NumLib::generateNodeValueDistribution(interpolate, *_msh, vec_node_ids);

    ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<float>::epsilon());
}


