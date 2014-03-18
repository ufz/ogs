/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <vector>
#include <memory>

#include "GEOObjects.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/LinearInterpolationAlongPolyline.h"
#include "MeshGeoToolsLib/LinearInterpolationAlongSurface.h"

#include "../TestTools.h"

class MeshGeoToolsLibInterpolationQuad : public testing::Test
{
public:
	MeshGeoToolsLibInterpolationQuad() :
		_geometric_size(10.0), _number_of_subdivisions_per_direction(10),
		_msh(MeshLib::MeshGenerator::generateRegularQuadMesh(_geometric_size, _number_of_subdivisions_per_direction)),
		_project_name("test"), _mshNodesSearcher(*_msh), _ply0(nullptr)
	{
		// create geometry
		std::vector<GeoLib::Point*>* pnts (new std::vector<GeoLib::Point*>);
		pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
		pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
		pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
		pnts->push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));

		std::vector<GeoLib::Polyline*>* plys (new std::vector<GeoLib::Polyline*>);
		_ply0 = new GeoLib::Polyline(*pnts);
		_ply0->addPoint(0);
		_ply0->addPoint(1);
		plys->push_back(_ply0);

		GeoLib::Polyline* ply1 = new GeoLib::Polyline(*pnts);
		ply1->addPoint(0);
		ply1->addPoint(1);
		ply1->addPoint(2);
		ply1->addPoint(3);
		ply1->addPoint(0);
		plys->push_back(ply1);

		std::vector<GeoLib::Surface*>* sfcs (new std::vector<GeoLib::Surface*>);
		_sfc1 = GeoLib::Surface::createSurface(*ply1);
		sfcs->push_back(_sfc1);

		_geo_objs.addPointVec(pnts,_project_name);
		_geo_objs.addPolylineVec(plys, _project_name);
		_geo_objs.addSurfaceVec(sfcs, _project_name);
	}

protected:
	const double _geometric_size;
	const std::size_t _number_of_subdivisions_per_direction;
	std::unique_ptr<MeshLib::Mesh> _msh;
	GeoLib::GEOObjects _geo_objs;
	std::string _project_name;
	MeshGeoTools::MeshNodeSearcher _mshNodesSearcher;
	GeoLib::Polyline* _ply0;
	GeoLib::Surface* _sfc1;
};

class MeshGeoToolsLibInterpolationHex : public testing::Test
{
public:
	MeshGeoToolsLibInterpolationHex() :
		_geometric_size(10.0), _number_of_subdivisions_per_direction(10),
		_msh(MeshLib::MeshGenerator::generateRegularHexMesh(_geometric_size, _number_of_subdivisions_per_direction)),
		_project_name("test"), _mshNodesSearcher(*_msh), _ply0(nullptr)
	{
		// create geometry
		std::vector<GeoLib::Point*>* pnts (new std::vector<GeoLib::Point*>);
		pnts->push_back(new GeoLib::Point(0.0, 0.0, 0.0));
		pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
		pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
		pnts->push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
		pnts->push_back(new GeoLib::Point(0.0, 0.0, _geometric_size));
		pnts->push_back(new GeoLib::Point(_geometric_size, 0.0, _geometric_size));
		pnts->push_back(new GeoLib::Point(_geometric_size, _geometric_size, _geometric_size));
		pnts->push_back(new GeoLib::Point(0.0, _geometric_size, _geometric_size));

		std::vector<GeoLib::Polyline*>* plys (new std::vector<GeoLib::Polyline*>);
		_ply0 = new GeoLib::Polyline(*pnts); // vertical polyline
		_ply0->addPoint(0);
		_ply0->addPoint(4);
		plys->push_back(_ply0);
		GeoLib::Polyline* ply1 = new GeoLib::Polyline(*pnts); // polygon for left surface
		ply1->addPoint(0);
		ply1->addPoint(3);
		ply1->addPoint(7);
		ply1->addPoint(4);
		ply1->addPoint(0);
		plys->push_back(ply1);

		std::vector<GeoLib::Surface*>* sfcs (new std::vector<GeoLib::Surface*>);
		_sfc1 = GeoLib::Surface::createSurface(*ply1);
		sfcs->push_back(_sfc1);

		_geo_objs.addPointVec(pnts,_project_name);
		_geo_objs.addPolylineVec(plys, _project_name);
		_geo_objs.addSurfaceVec(sfcs, _project_name);
	}

protected:
	const double _geometric_size;
	const std::size_t _number_of_subdivisions_per_direction;
	std::unique_ptr<MeshLib::Mesh> _msh;
	GeoLib::GEOObjects _geo_objs;
	std::string _project_name;
	MeshGeoTools::MeshNodeSearcher _mshNodesSearcher;
	GeoLib::Polyline* _ply0;
	GeoLib::Surface* _sfc1;
};

TEST_F(MeshGeoToolsLibInterpolationQuad, Polyline)
{
	const std::vector<std::size_t> vec_point_ids = {{0, 1}};
	const std::vector<double> vec_point_values = {{0., 100.}};
	std::vector<double> expected = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}};

	std::vector<double> interpolated_values;
	MeshGeoTools::LinearInterpolationAlongPolyline interpolate(
			_mshNodesSearcher.getMeshNodesAlongPolyline(*_ply0),
			vec_point_ids, vec_point_values,
			interpolated_values);

	ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(MeshGeoToolsLibInterpolationQuad, Surface)
{
	const std::vector<std::size_t> vec_point_ids = {{0, 1, 2, 3}};
	const std::vector<double> vec_point_values = {{0., 100., 100., 0.}};
	std::vector<double> expected(_msh->getNNodes());
	for (std::size_t i=0; i<_msh->getNNodes(); i++) {
		expected[i] = (i%(_number_of_subdivisions_per_direction+1)) * 10;
	}

	std::vector<double> interpolated_values;
	MeshGeoTools::LinearInterpolationAlongSurface interpolate(
			*_msh, _mshNodesSearcher.getMeshNodesAlongSurface(*_sfc1),
			vec_point_ids, vec_point_values, 0.0,
			interpolated_values);

	// the machine epsilon for double is too small for this test
	ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<float>::epsilon());
}

TEST_F(MeshGeoToolsLibInterpolationHex, Polyline)
{
	const std::vector<std::size_t> vec_point_ids = {{0, 4}};
	const std::vector<double> vec_point_values = {{0., 100.}};
	std::vector<double> expected = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}};

	std::vector<double> interpolated_values;
	MeshGeoTools::LinearInterpolationAlongPolyline interpolate(
			_mshNodesSearcher.getMeshNodesAlongPolyline(*_ply0),
			vec_point_ids, vec_point_values,
			interpolated_values);

	ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<double>::epsilon());
}

TEST_F(MeshGeoToolsLibInterpolationHex, Surface)
{
	const std::vector<std::size_t> vec_point_ids = {{0, 3, 7, 4}};
	const std::vector<double> vec_point_values = {{0., 100., 100., 0.}};
	std::vector<double> expected(std::pow(_number_of_subdivisions_per_direction+1, 2));
	for (std::size_t i=0; i<expected.size(); i++) {
		expected[i] = (i%(_number_of_subdivisions_per_direction+1)) * 10;
	}

	std::vector<double> interpolated_values;
	MeshGeoTools::LinearInterpolationAlongSurface interpolate(
			*_msh, _mshNodesSearcher.getMeshNodesAlongSurface(*_sfc1),
			vec_point_ids, vec_point_values, 0.0,
			interpolated_values);

	ASSERT_ARRAY_NEAR(expected, interpolated_values, expected.size(), std::numeric_limits<float>::epsilon());
}

