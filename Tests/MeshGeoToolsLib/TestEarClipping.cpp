/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "Applications/FileIO/Legacy/createSurface.h"
#include "GeoLib/EarClippingTriangulation.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/Polyline.h"
#include "InfoLib/TestInfo.h"

TEST(TestEarClipping, TestEarClippingPolygons)
{
    std::string const name =
        TestInfoLib::TestInfo::data_path + "/MeshGeoToolsLib/buildings.gml";
    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    ASSERT_TRUE(xml.readFile(name));
    auto geo_name = geo_objects.getGeometryNames()[0];
    auto polygons = *geo_objects.getPolylineVec(geo_name);
    ASSERT_EQ(4, polygons.size());
    for (auto polygon : polygons)
    {
        std::unique_ptr<GeoLib::Surface> sfc(
            FileIO::createSurfaceWithEarClipping(*polygon));
        ASSERT_TRUE(sfc != nullptr);
        ASSERT_EQ(sfc->getNumberOfTriangles(),
                  polygon->getNumberOfPoints() - 3);
    }
}
