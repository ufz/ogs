/**
 * \file   TestBoundingSphere.cpp
 * \author Karsten Rink
 * \date   2014-08-29
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "gtest/gtest.h"

#include "GeoLib/BoundingSphere.h"
#include "GeoLib/Point.h"

TEST(GeoLib, TestBoundingSphere)
{
    GeoLib::Point a(0,  0   ,0);
	GeoLib::Point b(2,  0   ,0);
	GeoLib::Point c(1,  0.1 ,0);
	GeoLib::Point d(1, -0.1 ,0);
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0,  0   , 0));
    pnts.push_back(new GeoLib::Point(2,  0   , 0));
    pnts.push_back(new GeoLib::Point(1,  0.1 , 0));
    pnts.push_back(new GeoLib::Point(1, -0.1 , 0));

    {
    GeoLib::BoundingSphere s(pnts);
	GeoLib::Point center = s.getCenter();
    ASSERT_NEAR(1.0, center[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.0, center[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.0, center[2], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(1.0, s.getRadius(), std::numeric_limits<double>::epsilon());
    }

    {
    (*pnts[2])[2] -= 1.4;
    GeoLib::BoundingSphere s(pnts);
	GeoLib::Point center = s.getCenter();
    ASSERT_NEAR(1.0, center[0], 0.0001);
    ASSERT_NEAR(0.0246, center[1], 0.0001);
    ASSERT_NEAR(-0.3446, center[2], 0.0001);
    ASSERT_NEAR(1.0580, s.getRadius(), 0.0001);
    }
    
    pnts[0] = new GeoLib::Point(0, 0, 0);
    pnts[1] = new GeoLib::Point(1, 0, 0);
    pnts[2] = new GeoLib::Point(1, 1, 0);
    pnts[3] = new GeoLib::Point(0, 1, 0);
    pnts.push_back(new GeoLib::Point(0, 0, 1));
    pnts.push_back(new GeoLib::Point(1, 0, 1));
    pnts.push_back(new GeoLib::Point(1, 1, 1));
    pnts.push_back(new GeoLib::Point(0, 1, 0.9));

    {
    GeoLib::BoundingSphere s(pnts);
	GeoLib::Point center = s.getCenter();
    ASSERT_NEAR(0.5, center[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.5, center[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.5, center[2], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.8660, s.getRadius(), 0.0001);
    }

    {
    (*pnts[7])[2] += 0.3;
    GeoLib::BoundingSphere s(pnts);
	GeoLib::Point center = s.getCenter();
    ASSERT_NEAR(0.5, center[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.5, center[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.6, center[2], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(0.9273, s.getRadius(), 0.0001);
    }
}
