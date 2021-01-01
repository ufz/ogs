/**
 * \date   2014-10-14
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>
#include <fstream>

#include "Applications/FileIO/Legacy/OGSIOVer4.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "InfoLib/TestInfo.h"
#include "filesystem.h"
#include "gtest/gtest.h"

class OGSIOVer4InterfaceTest : public ::testing::Test
{
public:
    OGSIOVer4InterfaceTest()
        : _test_path(fs::temp_directory_path() /= BaseLib::randomString(32)),
          _gli_fname(_test_path), _surface_fname(_test_path)
    {
        fs::create_directory(_test_path);
        std::ofstream gli_out(_gli_fname /= "test.gli");
        gli_out << "#POINTS\n";
        gli_out << "0 0 0 0\n";
        gli_out << "#SURFACE\n";
        gli_out << "  $NAME\n";
        gli_out << "    Surface\n";
        gli_out << "  $TIN\n";
        gli_out << "    Surface.tin\n";
        gli_out << "#STOP\n";
        gli_out.close();
        _surface_fname /= "Surface.tin";
    }

    ~OGSIOVer4InterfaceTest() override { fs::remove_all(_test_path); }

protected:
    const fs::path _test_path;
    fs::path _gli_fname;
    fs::path _surface_fname;
};

TEST_F(OGSIOVer4InterfaceTest, SimpleTIN)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0 0.0 1.0 1.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs != nullptr);
    ASSERT_EQ(1u, geometries.getSurfaceVec(geometry_name)->size());
    ASSERT_EQ(2u, (*geometries.getSurfaceVec(geometry_name))[0]->getNumberOfTriangles());
}

TEST_F(OGSIOVer4InterfaceTest, StillCorrectTINWihtAdditionalValueAtEndOfLine)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0 10\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs != nullptr);
    ASSERT_EQ(1u, geometries.getSurfaceVec(geometry_name)->size());
    ASSERT_EQ(2u, (*geometries.getSurfaceVec(geometry_name))[0]->getNumberOfTriangles());
}

TEST_F(OGSIOVer4InterfaceTest, InvalidTIN_ZeroAreaTri)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 0.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs == nullptr);
}


TEST_F(OGSIOVer4InterfaceTest, InvalidTIN_LineDoesNotStartWithID)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
    tin_out << "a\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs == nullptr);
}


TEST_F(OGSIOVer4InterfaceTest, InvalidTIN_PointIsMissing)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs == nullptr);
}

TEST_F(OGSIOVer4InterfaceTest, InvalidTIN_CoordOfPointIsMissing)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs == nullptr);
}

TEST_F(OGSIOVer4InterfaceTest, SimpleTIN_AdditionalEmptyLinesAtEnd)
{
    std::ofstream tin_out (_surface_fname);
    tin_out << "0 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0 10\n";
    tin_out << "1 0.0 0.0 0.0 1.0 0.0.0 0.0 0.0 1.0\n\n\n";
    tin_out.close();

    // read geometry
    GeoLib::GEOObjects geometries;
    std::vector<std::string> errors;
    std::string geometry_name("TestGeometry");
    FileIO::Legacy::readGLIFileV4(_gli_fname.string(), geometries, geometry_name,
                                  errors, "dummy_for_gmsh_path");

    std::vector<GeoLib::Surface*> const*
        sfcs(geometries.getSurfaceVec(geometry_name));
    ASSERT_TRUE(sfcs != nullptr);
    ASSERT_EQ(1u, geometries.getSurfaceVec(geometry_name)->size());
    ASSERT_EQ(2u, (*geometries.getSurfaceVec(geometry_name))[0]->getNumberOfTriangles());
}
