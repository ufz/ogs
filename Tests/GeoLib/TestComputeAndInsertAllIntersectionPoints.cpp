/**
 * @file TestComputeAndInsertAllIntersectionPoints.cpp
 * @date 2014-05-09
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <ctime>
#include <tuple>

#include "gtest/gtest.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"

TEST(GeoLib, TestComputeAndInsertAllIntersectionPoints)
{
    GeoLib::GEOObjects geo_objs;
    std::string geo_name("TestGeometry");

    {
        // *** insert points in vector
        auto pnts = std::make_unique<std::vector<GeoLib::Point*>>();
        pnts->push_back(new GeoLib::Point(0.0,0.0,0.0,0));
        pnts->push_back(new GeoLib::Point(11.0,0.0,0.0,1));

        pnts->push_back(new GeoLib::Point(0.0,1.0,0.0,2));
        pnts->push_back(new GeoLib::Point(1.0,-1.0,0.0,3));
        pnts->push_back(new GeoLib::Point(2.0, 1.0,0.0,4));
        pnts->push_back(new GeoLib::Point(3.0,-1.0,0.0,5));
        pnts->push_back(new GeoLib::Point(4.0, 1.0,0.0,6));
        pnts->push_back(new GeoLib::Point(5.0,-1.0,0.0,7));
        pnts->push_back(new GeoLib::Point(6.0, 1.0,0.0,8));
        pnts->push_back(new GeoLib::Point(7.0,-1.0,0.0,9));
        pnts->push_back(new GeoLib::Point(8.0, 1.0,0.0,10));
        pnts->push_back(new GeoLib::Point(9.0,-1.0,0.0,11));
        pnts->push_back(new GeoLib::Point(10.0, 1.0,0.0,12));
        pnts->push_back(new GeoLib::Point(11.0,-1.0,0.0,13));

        geo_objs.addPointVec(std::move(pnts), geo_name);
    }

    // *** create polylines
    auto& pnts = *geo_objs.getPointVec(geo_name);
    auto* ply0(new GeoLib::Polyline(pnts));
    ply0->addPoint(0);
    ply0->addPoint(1);
    auto* ply1(new GeoLib::Polyline(pnts));
    for (std::size_t k(2); k<pnts.size(); ++k)
        ply1->addPoint(k);
    auto* plys(new std::vector<GeoLib::Polyline*>);
    plys->push_back(ply0);
    plys->push_back(ply1);

    GeoLib::PointVec &pnt_vec(*(const_cast<GeoLib::PointVec*>(geo_objs.getPointVecObj(geo_name))));
    GeoLib::computeAndInsertAllIntersectionPoints(pnt_vec, *plys);

    ASSERT_EQ(25u, pnt_vec.size());
    ASSERT_EQ(13u, ply0->getNumberOfPoints());
    ASSERT_EQ(23u, ply1->getNumberOfPoints());

    // check correct order of points in ply0
    EXPECT_EQ(0u, ply0->getPointID(0));
    EXPECT_EQ(14u, ply0->getPointID(1));
    EXPECT_EQ(15u, ply0->getPointID(2));
    EXPECT_EQ(16u, ply0->getPointID(3));
    EXPECT_EQ(17u, ply0->getPointID(4));
    EXPECT_EQ(18u, ply0->getPointID(5));
    EXPECT_EQ(19u, ply0->getPointID(6));
    EXPECT_EQ(20u, ply0->getPointID(7));
    EXPECT_EQ(21u, ply0->getPointID(8));
    EXPECT_EQ(22u, ply0->getPointID(9));
    EXPECT_EQ(23u, ply0->getPointID(10));
    EXPECT_EQ(24u, ply0->getPointID(11));
    EXPECT_EQ(1u, ply0->getPointID(12));

    // check correct order of points in ply1
    EXPECT_EQ(2u, ply1->getPointID(0));
    EXPECT_EQ(14u, ply1->getPointID(1));
    EXPECT_EQ(3u, ply1->getPointID(2));
    EXPECT_EQ(15u, ply1->getPointID(3));
    EXPECT_EQ(4u, ply1->getPointID(4));
    EXPECT_EQ(16u, ply1->getPointID(5));
    EXPECT_EQ(5u, ply1->getPointID(6));
    EXPECT_EQ(17u, ply1->getPointID(7));
    EXPECT_EQ(6u, ply1->getPointID(8));
    EXPECT_EQ(18u, ply1->getPointID(9));
    EXPECT_EQ(7u, ply1->getPointID(10));
    EXPECT_EQ(19u, ply1->getPointID(11));
    EXPECT_EQ(8u, ply1->getPointID(12));
    EXPECT_EQ(20u, ply1->getPointID(13));
    EXPECT_EQ(9u, ply1->getPointID(14));
    EXPECT_EQ(21u, ply1->getPointID(15));
    EXPECT_EQ(10u, ply1->getPointID(16));
    EXPECT_EQ(22u, ply1->getPointID(17));
    EXPECT_EQ(11u, ply1->getPointID(18));
    EXPECT_EQ(23u, ply1->getPointID(19));
    EXPECT_EQ(12u, ply1->getPointID(20));
    EXPECT_EQ(24u, ply1->getPointID(21));
    EXPECT_EQ(13u, ply1->getPointID(22));

    delete plys;
    delete ply1;
    delete ply0;
}

