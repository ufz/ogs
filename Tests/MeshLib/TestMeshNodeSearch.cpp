/**
 * @file TestMeshNodeSearch.cpp
 * @author git blame TestMeshNodeSearch.cpp
 * @date Oct 28, 2013
 * @brief Test the implementation of class MeshNodeSearch.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshGenerators/MeshGenerator.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

using namespace MeshLib;

class MeshLibMeshNodeSearchInSimpleQuadMesh : public testing::Test
{
public:
	MeshLibMeshNodeSearchInSimpleQuadMesh() :
		_geometric_size(10.0), _number_of_subdivisions_per_direction(99),
		_quad_mesh(MeshGenerator::generateRegularQuadMesh(_geometric_size, _number_of_subdivisions_per_direction))
	{}

	~MeshLibMeshNodeSearchInSimpleQuadMesh()
	{
		delete _quad_mesh;
	}

protected:
	const double _geometric_size;
	const std::size_t _number_of_subdivisions_per_direction;
	Mesh* _quad_mesh;
};

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PointSearch)
{
	ASSERT_TRUE(_quad_mesh != nullptr);
	// 1 create a geometry
	GeoLib::Point pnt(0.0, 0.0, 0.0);

	// 2 perform search and compare results with expected vals
	MeshGeoTools::MeshNodeSearcher mesh_node_searcher(*_quad_mesh);

	// find ORIGIN
	ASSERT_EQ(0, mesh_node_searcher.getMeshNodeIDForPoint(pnt));

	pnt[0] = 0.049;
	pnt[1] = 0.049;
	ASSERT_EQ(0, mesh_node_searcher.getMeshNodeIDForPoint(pnt));

	pnt[0] = 0.051;
	pnt[1] = 0.049;
	ASSERT_EQ(1, mesh_node_searcher.getMeshNodeIDForPoint(pnt));

	pnt[0] = 0.049;
	pnt[1] = 0.051;
	ASSERT_EQ(100, mesh_node_searcher.getMeshNodeIDForPoint(pnt));

	pnt[0] = 0.051;
	pnt[1] = 0.051;
	ASSERT_EQ(101, mesh_node_searcher.getMeshNodeIDForPoint(pnt));

	pnt[0] = 9.951;
	pnt[1] = 9.951;
	ASSERT_EQ((_number_of_subdivisions_per_direction+1) * (_number_of_subdivisions_per_direction+1) - 1, mesh_node_searcher.getMeshNodeIDForPoint(pnt));
}

TEST_F(MeshLibMeshNodeSearchInSimpleQuadMesh, PolylineSearch)
{
	ASSERT_TRUE(_quad_mesh != nullptr);
	// create geometry
	std::vector<GeoLib::Point*> pnts;
	pnts.push_back(new GeoLib::Point(0.0, 0.0, 0.0));
	pnts.push_back(new GeoLib::Point(_geometric_size, 0.0, 0.0));
	pnts.push_back(new GeoLib::Point(0.0, 0.049, 0.0));
	pnts.push_back(new GeoLib::Point(_geometric_size, 0.049, 0.0));
	pnts.push_back(new GeoLib::Point(0.0, _geometric_size, 0.0));
	pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size, 0.0));
	pnts.push_back(new GeoLib::Point(0.0, _geometric_size - 0.049, 0.0));
	pnts.push_back(new GeoLib::Point(_geometric_size, _geometric_size - 0.049, 0.0));

	GeoLib::Polyline ply0(pnts);
	ply0.addPoint(0);
	ply0.addPoint(1);

	// perform search and compare results with expected vals
	MeshGeoTools::MeshNodeSearcher mesh_node_searcher(*_quad_mesh);
	std::vector<std::size_t> const& found_ids_ply0(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply0));

	ASSERT_EQ(100, found_ids_ply0.size());
	for (std::size_t k(0); k<found_ids_ply0.size(); k++)
		ASSERT_EQ(k, found_ids_ply0[k]);

	GeoLib::Polyline ply1(pnts);
	ply1.addPoint(2);
	ply1.addPoint(3);
	std::vector<std::size_t> const& found_ids_ply1(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply1));

	ASSERT_EQ(100, found_ids_ply1.size());
	for (std::size_t k(0); k<found_ids_ply1.size(); k++)
		ASSERT_EQ(k, found_ids_ply1[k]);

	GeoLib::Polyline ply2(pnts);
	ply2.addPoint(4);
	ply2.addPoint(5);
	std::vector<std::size_t> const& found_ids_ply2(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply2));

	std::size_t offset((_number_of_subdivisions_per_direction+1)*_number_of_subdivisions_per_direction);
	ASSERT_EQ(100, found_ids_ply2.size());
	for (std::size_t k(0); k<found_ids_ply2.size(); k++)
		ASSERT_EQ(offset + k, found_ids_ply2[k]);

	GeoLib::Polyline ply3(pnts);
	ply3.addPoint(6);
	ply3.addPoint(7);
	std::vector<std::size_t> const& found_ids_ply3(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply3));

	ASSERT_EQ(100, found_ids_ply3.size());
	for (std::size_t k(0); k<found_ids_ply3.size(); k++)
		ASSERT_EQ(offset + k, found_ids_ply3[k]);

	// left border
	GeoLib::Polyline ply4(pnts);
	ply4.addPoint(0);
	ply4.addPoint(6);
	std::vector<std::size_t> const& found_ids_ply4(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply4));

	ASSERT_EQ(100, found_ids_ply4.size());
	for (std::size_t k(0); k<found_ids_ply4.size(); k++)
		ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1), found_ids_ply4[k]);

	// right border
	GeoLib::Polyline ply5(pnts);
	ply5.addPoint(1);
	ply5.addPoint(7);
	std::vector<std::size_t> const& found_ids_ply5(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply5));

	ASSERT_EQ(100, found_ids_ply5.size());
	for (std::size_t k(0); k<found_ids_ply5.size(); k++)
		ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1)+_number_of_subdivisions_per_direction, found_ids_ply5[k]);

	// diagonal
	GeoLib::Polyline ply6(pnts);
	ply6.addPoint(0);
	ply6.addPoint(5);
	std::vector<std::size_t> const& found_ids_ply6(mesh_node_searcher.getMeshNodeIDsAlongPolyline(ply6));

	ASSERT_EQ(100, found_ids_ply6.size());
	for (std::size_t k(0); k<found_ids_ply6.size(); k++)
		ASSERT_EQ(k*(_number_of_subdivisions_per_direction+1)+k, found_ids_ply6[k]);

	std::for_each(pnts.begin(), pnts.end(), [](GeoLib::Point* pnt) { delete pnt; });
}
