// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ctime>

#include "MathLib/Point3dWithID.h"

using namespace MathLib;

TEST(MathLib, Point3dWithID)
{
    Point3dWithID p0(0, 0, 0, 1);
    const Point3dWithID& p1(p0);  // copy constructor
    Point3dWithID p2(p0, 2);      // constructor for resetting the id

    EXPECT_EQ(p0.getID(), p1.getID());
    EXPECT_NE(p0.getID(), p2.getID());
}
