/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "FileIO/XmlIO/Qt/XmlPrjInterface.h"

TEST(TestPrjInterface, QtPrjWriterReaderTest)
{
	std::string const orig_project_file(BaseLib::BuildInfo::data_path + 
			"/Elliptic/cube_1x1x1_GroundWaterFlow/cube_1e4.prj");
	std::string test_project_file = BaseLib::BuildInfo::tests_tmp_path + "TestProjectIO.prj";

	std::unique_ptr<MeshLib::Mesh> comp_mesh;
	std::vector<GeoLib::Point*> comp_points;

	{
		// read project from test repo
		ProjectData project;
		FileIO::XmlPrjInterface xml(project);
		xml.readFile(orig_project_file);
		std::vector<MeshLib::Mesh*> const& mesh_vec = project.getMeshObjects();
		ASSERT_EQ(1, mesh_vec.size());
		ASSERT_EQ(12167, mesh_vec[0]->getNNodes());
		ASSERT_EQ(10648, mesh_vec[0]->getNElements());
		std::vector<std::string> geo_names;
		project.getGEOObjects()->getGeometryNames(geo_names);
		ASSERT_EQ(1, geo_names.size());
		std::vector<GeoLib::Point*> const& points (*project.getGEOObjects()->getPointVec(geo_names[0]));
		ASSERT_EQ(8, points.size());

		// copy data for comparison later
		comp_mesh.reset(new MeshLib::Mesh(*mesh_vec[0]));
		for (std::size_t i=0; i<points.size(); ++i)
			comp_points.push_back(new GeoLib::Point(*points[i]));

		// write project to temp
		ASSERT_EQ(0, xml.writeToFile(test_project_file));
	}

	// read project from temp again
	ProjectData project;
	FileIO::XmlPrjInterface xml(project);
	xml.readFile(test_project_file);
	std::vector<MeshLib::Mesh*> const& mesh_vec = project.getMeshObjects();
	ASSERT_EQ(1, mesh_vec.size());
	ASSERT_EQ(comp_mesh->getNNodes(), mesh_vec[0]->getNNodes());
	ASSERT_EQ(comp_mesh->getNElements(), mesh_vec[0]->getNElements());
	std::vector<std::string> geo_names;
	project.getGEOObjects()->getGeometryNames(geo_names);
	ASSERT_EQ(1, geo_names.size());
	std::vector<GeoLib::Point*> const& points (*project.getGEOObjects()->getPointVec(geo_names[0]));
	ASSERT_EQ(8, points.size());
	for (std::size_t i=0; i<points.size(); ++i)
	{
		ASSERT_EQ((*comp_points[i])[0], (*points[i])[0]);
		ASSERT_EQ((*comp_points[i])[1], (*points[i])[1]);
		ASSERT_EQ((*comp_points[i])[2], (*points[i])[2]);
	}

	// remove data from temp
	std::remove(test_project_file.c_str());
	test_project_file += ".md5";
	std::remove(test_project_file.c_str());
	std::string geo_name = geo_names[0] + ".gml";
	std::remove(geo_name.c_str());
	geo_name += ".md5";
	std::remove(geo_name.c_str());
	std::string mesh_name (mesh_vec[0]->getName() + ".vtu");
	std::remove(mesh_name.c_str());
}
