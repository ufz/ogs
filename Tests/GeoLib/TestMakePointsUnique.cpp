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

// Testing nullptr input vector.
TEST(GeoLib, TestPointVecCtorNullptr)
{
	ASSERT_ANY_THROW(GeoLib::PointVec("JustAName", nullptr));
}

// Testing empty input vector.
TEST(GeoLib, TestPointVecCtorEmpty)
{
	std::vector<GeoLib::Point*> ps;
	ASSERT_ANY_THROW(GeoLib::PointVec("JustAName", &ps));
}

// Testing input vector with single point.
TEST(GeoLib, TestPointVecCtorSinglePoint)
{
	std::vector<GeoLib::Point*> ps;
	ps.push_back(new GeoLib::Point(0,0,0));
	FAIL() << "SEGV Error in destructor TemplateVec.h:62.\n"
		<< "Dealloction of the points vector fails.";
	ASSERT_NO_THROW(GeoLib::PointVec point_vec("JustAName", &ps));
}

}
