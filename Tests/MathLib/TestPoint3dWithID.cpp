/**
 * \file
 * \date 2015-05-21
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <ctime>

#include "MathLib/Point3dWithID.h"
#include "gtest/gtest.h"

using namespace MathLib;

TEST(MathLib, Point3dWithID)
{
    Point3dWithID p0(0, 0, 0, 1);
    const Point3dWithID& p1(p0);  // copy constructor
    Point3dWithID p2(p0, 2);      // constructor for resetting the id

    EXPECT_EQ(p0.getID(), p1.getID());
    EXPECT_NE(p0.getID(), p2.getID());
}
