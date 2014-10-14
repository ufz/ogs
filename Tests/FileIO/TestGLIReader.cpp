/**
 * \date   2014-10-14
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"

// FileIO
#include "FileIO/Legacy/OGSIOVer4.h"

// GeoLib
#include "GeoLib/GEOObjects.h"

class OGSIOVer4InterfaceTest : public ::testing::Test
{
public:
	OGSIOVer4InterfaceTest()
		: _gli_fname(BaseLib::BuildInfo::tests_tmp_path+"test.gli")
	{
		std::ofstream gli_out(_gli_fname);
		gli_out << "#POINTS\n";
		gli_out << "0 0 0 0\n";
		gli_out << "#SURFACE\n";
		gli_out << "  $NAME\n";
		gli_out << "    Surface\n";
		gli_out << "  $TIN\n";
		gli_out << "    Surface.tin\n";
		gli_out << "#STOP\n";
		gli_out.close();
	}

	~OGSIOVer4InterfaceTest()
	{
		boost::filesystem::remove(_gli_fname);
	}

protected:
	std::string _gli_fname;
};

TEST_F(OGSIOVer4InterfaceTest, SimpleTIN)
{
	std::string tin_fname(BaseLib::BuildInfo::tests_tmp_path+"Surface.tin");
	std::ofstream tin_out (tin_fname);
	tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
	tin_out << "1 0.0 0.0 0.0 1.0 0.0.0 0.0 1.0 1.0\n";
	tin_out.close();

	// read geometry
	GeoLib::GEOObjects geometries;
	std::vector<std::string> errors;
	std::string geometry_name("TestGeometry");
	FileIO::Legacy::readGLIFileV4(_gli_fname, &geometries, geometry_name, errors);

	std::vector<GeoLib::Surface*> const*
		sfcs(geometries.getSurfaceVec(geometry_name));
	ASSERT_TRUE(sfcs != nullptr);
	ASSERT_EQ(1u, geometries.getSurfaceVec(geometry_name)->size());
	ASSERT_EQ(2u, (*geometries.getSurfaceVec(geometry_name))[0]->getNTriangles());

	boost::filesystem::remove(tin_fname);
}

TEST_F(OGSIOVer4InterfaceTest, InvalidTIN)
{
	std::string tin_fname(BaseLib::BuildInfo::tests_tmp_path+"Surface.tin");
	std::ofstream tin_out (tin_fname);
	tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
	tin_out << "1 0.0 0.0 0.0 1.0 0.0.0\n";
	tin_out.close();

	// read geometry
	GeoLib::GEOObjects geometries;
	std::vector<std::string> errors;
	std::string geometry_name("TestGeometry");
	FileIO::Legacy::readGLIFileV4(_gli_fname, &geometries, geometry_name, errors);

	std::vector<GeoLib::Surface*> const*
		sfcs(geometries.getSurfaceVec(geometry_name));
	ASSERT_TRUE(sfcs == nullptr);

	boost::filesystem::remove(tin_fname);
}

TEST_F(OGSIOVer4InterfaceTest, InvalidTIN_II)
{
	std::string tin_fname(BaseLib::BuildInfo::tests_tmp_path+"Surface.tin");
	std::ofstream tin_out (tin_fname);
	tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0\n";
	tin_out << "1 0.0 0.0 0.0 1.0 0.0.0\n";
	tin_out.close();

	// read geometry
	GeoLib::GEOObjects geometries;
	std::vector<std::string> errors;
	std::string geometry_name("TestGeometry");
	FileIO::Legacy::readGLIFileV4(_gli_fname, &geometries, geometry_name, errors);

	std::vector<GeoLib::Surface*> const*
		sfcs(geometries.getSurfaceVec(geometry_name));
	ASSERT_TRUE(sfcs == nullptr);

	boost::filesystem::remove(tin_fname);
}

