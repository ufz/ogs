/**
 * \file   TestXmlGmlReader.cpp
 * \author Karsten Rink
 * \date   2013-03-20
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"
#include "Configure.h"

// FileIO
#include "XmlIO/Qt/XmlGmlInterface.h"

// OgsLib
#include "OGS/ProjectData.h"

// GeoLib
#include "Polyline.h"


TEST(FileIO, XmlGmlWriterReaderTest)
{
	// Writer test
	std::string test_data_file(std::string(SOURCEPATH).append("/Tests/FileIO/xmlgmltestdata.gml"));

	ProjectData project;
	GeoLib::GEOObjects geo_objects;

	//setup test data
	std::string geo_name("TestData");
	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>(10);
	std::vector<GeoLib::Polyline*> *lines = new std::vector<GeoLib::Polyline*>(5);
	std::vector<GeoLib::Surface*> *sfcs = new std::vector<GeoLib::Surface*>(2);
	std::map<std::string, std::size_t>* ply_names = new std::map<std::string, std::size_t>;

	(*points)[0] = new GeoLib::Point(1,1,0);
	(*points)[1] = new GeoLib::Point(1,1,0);
	(*points)[2] = new GeoLib::Point(1,2,0);
	(*points)[3] = new GeoLib::Point(1,3,0);
	(*points)[4] = new GeoLib::Point(2,1,0);
	(*points)[5] = new GeoLib::Point(2,2,0);
	(*points)[6] = new GeoLib::Point(2,3,0);
	(*points)[7] = new GeoLib::Point(3,1,0);
	(*points)[8] = new GeoLib::Point(3,2,0);
	(*points)[9] = new GeoLib::Point(3,3,0);
	geo_objects.addPointVec(points, geo_name);
	const std::vector<std::size_t> pnt_id_map (geo_objects.getPointVecObj(geo_name)->getIDMap());

	(*lines)[0] = new GeoLib::Polyline(*points);
	(*lines)[0]->addPoint(pnt_id_map[0]); (*lines)[0]->addPoint(pnt_id_map[2]); (*lines)[0]->addPoint(pnt_id_map[3]);
	ply_names->insert(std::pair<std::string, std::size_t>("left", 0));
	(*lines)[1] = new GeoLib::Polyline(*points);
	(*lines)[1]->addPoint(pnt_id_map[4]); (*lines)[1]->addPoint(pnt_id_map[5]); (*lines)[1]->addPoint(pnt_id_map[6]);
	ply_names->insert(std::pair<std::string, std::size_t>("center", 1));
	(*lines)[2] = new GeoLib::Polyline(*points);
	(*lines)[2]->addPoint(pnt_id_map[1]); (*lines)[2]->addPoint(pnt_id_map[4]);
	(*lines)[3] = new GeoLib::Polyline(*points);
	(*lines)[3]->addPoint(pnt_id_map[4]); (*lines)[3]->addPoint(pnt_id_map[7]);
	(*lines)[4] = new GeoLib::Polyline(*points);
	(*lines)[4]->addPoint(pnt_id_map[7]); (*lines)[4]->addPoint(pnt_id_map[8]); (*lines)[4]->addPoint(pnt_id_map[9]);
	ply_names->insert(std::pair<std::string, std::size_t>("right", 4));
	geo_objects.addPolylineVec(lines, geo_name, ply_names);

	(*sfcs)[0] = new GeoLib::Surface(*points);
	(*sfcs)[0]->addTriangle(pnt_id_map[1],pnt_id_map[4],pnt_id_map[2]);
	(*sfcs)[0]->addTriangle(pnt_id_map[2],pnt_id_map[4],pnt_id_map[5]);
	(*sfcs)[0]->addTriangle(pnt_id_map[2],pnt_id_map[5],pnt_id_map[3]);
	(*sfcs)[0]->addTriangle(pnt_id_map[3],pnt_id_map[5],pnt_id_map[6]);
	(*sfcs)[1] = new GeoLib::Surface(*points);
	(*sfcs)[1]->addTriangle(pnt_id_map[4],pnt_id_map[7],pnt_id_map[9]);
	(*sfcs)[1]->addTriangle(pnt_id_map[4],pnt_id_map[9],pnt_id_map[6]);
	geo_objects.addSurfaceVec(sfcs, geo_name);

	FileIO::XmlGmlInterface xml(geo_objects);
	xml.setNameForExport(geo_name);
	int result = xml.writeToFile(test_data_file);
	ASSERT_EQ(result, 1);

	// Reader test
	result = xml.readFile(QString::fromStdString(test_data_file));
	ASSERT_EQ(result, 1);

	const std::vector<GeoLib::Point*> *readerPoints = geo_objects.getPointVec(geo_name);
	const GeoLib::PolylineVec *line_vec = geo_objects.getPolylineVecObj(geo_name);
	const std::vector<GeoLib::Polyline*> *readerLines = geo_objects.getPolylineVec(geo_name);
	const std::vector<GeoLib::Surface*> *readerSfcs = geo_objects.getSurfaceVec(geo_name);
	ASSERT_EQ(9u, readerPoints->size());
	ASSERT_EQ(5u, readerLines->size());
	ASSERT_EQ(2u, readerSfcs->size());

	GeoLib::Point* pnt = (*readerPoints)[7];
	ASSERT_EQ(3.0, (*pnt)[0]);
	ASSERT_EQ(2.0, (*pnt)[1]);
	ASSERT_EQ(0.0, (*pnt)[2]);

	GeoLib::Polyline* line = (*readerLines)[4];
	ASSERT_EQ(3u, line->getNumberOfPoints());
	ASSERT_EQ(6u, line->getPointID(0));
	ASSERT_EQ(7u, line->getPointID(1));
	ASSERT_EQ(8u, line->getPointID(2));
	std::string line_name("");
	line_vec->getNameOfElementByID(4, line_name);
	ASSERT_EQ("right", line_name);

	GeoLib::Surface* sfc = (*readerSfcs)[1];
	ASSERT_EQ(2u, sfc->getNTriangles());
	const GeoLib::Triangle* tri = (*sfc)[1];
	ASSERT_EQ(3u, (*tri)[0]);
	ASSERT_EQ(8u, (*tri)[1]);
	ASSERT_EQ(5u, (*tri)[2]);

	boost::filesystem::remove(test_data_file);
	test_data_file += ".md5";
	boost::filesystem::remove(test_data_file);
}
