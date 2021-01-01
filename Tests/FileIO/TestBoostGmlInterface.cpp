/**
 * \file
 * \author Karsten Rink
 * \date   2016-02-16
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

#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/TestInfo.h"
#include "Tests/FileIO/TestGmlInterface.h"
#include "filesystem.h"
#include "gtest/gtest.h"

TEST_F(TestGmlInterface, BoostXmlGmlWriterReaderTest)
{
    // Writer test
    std::string test_data_file =
        (fs::temp_directory_path() /= BaseLib::randomString(32) + ".gml").string();

    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    xml.setNameForExport(geo_name);
    int result = xml.writeToFile(test_data_file);
    EXPECT_EQ(result, 1);

    // remove the written data from the data structures
    geo_objects.removeSurfaceVec(geo_name);
    geo_objects.removePolylineVec(geo_name);
    geo_objects.removePointVec(geo_name);

    // Reader test
    result = xml.readFile(test_data_file);
    EXPECT_EQ(1, result);

    std::remove(test_data_file.c_str());
    test_data_file += ".md5";
    std::remove(test_data_file.c_str());

    checkPointProperties();
    checkPolylineProperties();
    checkSurfaceProperties();
}
