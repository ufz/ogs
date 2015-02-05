/**
 * @file TestMeshValidation.cpp
 * @author Karsten Rink
 * @date 2015-01-29
 * @brief Tests for MeshRevision class
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/ConvertRasterToMesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/Elements/Element.h"

#include "GeoLib/Raster.h"

TEST(MeshValidation, UnusedNodes)
{
	std::array<double, 12> pix = {0,0.1,0.2,0.1,0,0,0.1,0,0,0,-0.1,0};
	GeoLib::Raster raster(4,3,0,0,1,pix.begin(), pix.end());
	MeshLib::ConvertRasterToMesh conv(raster, MeshElemType::TRIANGLE, MeshLib::UseIntensityAs::ELEVATION);
	MeshLib::Mesh* mesh = conv.execute();
	std::vector<std::size_t> u_nodes = MeshLib::MeshValidation::findUnusedMeshNodes(*mesh);
	ASSERT_EQ(0, u_nodes.size());

	std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
	nodes.push_back(new MeshLib::Node(-1,-1,-1));
	std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
	MeshLib::Mesh mesh2("mesh2", nodes, elems);
	u_nodes = MeshLib::MeshValidation::findUnusedMeshNodes(mesh2);
	ASSERT_EQ(1, u_nodes.size());
	ASSERT_EQ(nodes.back()->getID(), u_nodes[0]);

	delete mesh;
}

TEST(MeshValidation, DetectHolesTri)
{
	std::array<double, 12> pix = {0,0.1,0.2,0.1,0,0,0.1,0,0,0,-0.1,0};
	GeoLib::Raster raster(4,3,0,0,1,pix.begin(), pix.end());
	MeshLib::ConvertRasterToMesh conv(raster, MeshElemType::TRIANGLE, MeshLib::UseIntensityAs::ELEVATION);
	MeshLib::Mesh* mesh = conv.execute();
	unsigned n_holes = MeshLib::MeshValidation::detectHoles(*mesh);
	ASSERT_EQ(0, n_holes);

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+12);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(1, n_holes);
	}

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+11);
		elems.erase(elems.begin()+11);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(1, n_holes);
	}

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+10);
		elems.erase(elems.begin()+12);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(2, n_holes);
	}

	delete mesh;
}

TEST(MeshValidation, DetectHolesHex)
{
	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(5,4,4,1,1,1,GeoLib::ORIGIN, "mesh");
	unsigned n_holes = MeshLib::MeshValidation::detectHoles(*mesh);
	ASSERT_EQ(0, n_holes);

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+27);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(1, n_holes);
	}

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+28);
		elems.erase(elems.begin()+27);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(1, n_holes);
	}

	{
		std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
		std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
		elems.erase(elems.begin()+29);
		elems.erase(elems.begin()+27);
		MeshLib::Mesh mesh2("mesh2", nodes, elems);
		n_holes = MeshLib::MeshValidation::detectHoles(mesh2);
		ASSERT_EQ(1, n_holes);
	}
	delete mesh;
}




