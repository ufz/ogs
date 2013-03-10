/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "GeoLib/PointVec.h"


// Testing empty input vector.
TEST(GeoLib, TestPointVecCtorEmpty)
{
	std::vector<GeoLib::Point*> ps;
	ASSERT_ANY_THROW(GeoLib::PointVec("JustAName", &ps));
}
