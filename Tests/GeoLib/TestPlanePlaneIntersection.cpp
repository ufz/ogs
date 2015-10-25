/**
 * @date 2015-10-19
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <limits>

#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>

#include "Tests/MathLib/AutoCheckTools.h"

#include "MathLib/MathTools.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/Plane.h"

namespace ac = autocheck;

struct GeoLibComputePlanePlaneIntersection : public ::testing::Test
{
	// to generate arbitrary straight line segments
	ac::randomTupleGenerator<double, 3> tuple_generator;
	ac::cons_generator<MathLib::Point3d, ac::randomTupleGenerator<double, 3>>
	    points_gen{tuple_generator};

	ac::gtest_reporter gtest_reporter;

	// First spanning vector of the planes must be coincident with d0.
	bool check(GeoLib::Plane plane1,
	           GeoLib::Plane plane2,
	           MathLib::Vector3 const& d0)
	{
		// reconstructed intersection line
		MathLib::Vector3 d;
		MathLib::Point3d p;

		std::tie(d, p) = GeoLib::computePlanePlaneIntersection(plane1, plane2);

		// check if the given vector d0 and the computed vector are parallel
		if (GeoLib::isParallel(d, d0) && plane1.isPointInPlane(p)
			&& plane2.isPointInPlane(p))
			return true;

		return false;
	}
};

TEST_F(GeoLibComputePlanePlaneIntersection, TestPlanePlaneIntersection)
{
	auto checkPlanePlaneIntersection =
		[this](std::vector<MathLib::Point3d> & pnts) -> bool
	{
		// given intersection line
		MathLib::Vector3 const d0(pnts[0]);
		MathLib::Vector3 const p0(pnts[1]);

		// Both planes go through p0
		return this->check({d0, pnts[2], p0}, {d0, pnts[3], p0}, d0);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(4,list_of(points_gen))),
		gtest_reporter);
}

TEST_F(GeoLibComputePlanePlaneIntersection,
	TestPlaneVerticalPlaneIntersection)
{
	auto checkPlanePlaneIntersection =
		[this](std::vector<MathLib::Point3d> & pnts) -> bool
	{
		// given intersection line
		MathLib::Vector3 const d0(pnts[0]);
		MathLib::Vector3 const p0(pnts[1]);

		// Both planes go through p0; First plane is vertical.
		return this->check({d0, {0.0, 0.0, 1.0}, p0}, {d0, pnts[2], p0}, d0);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(3,list_of(points_gen))),
		gtest_reporter);
}

TEST_F(GeoLibComputePlanePlaneIntersection,
	TestHorizontalPlaneVerticalPlaneIntersection)
{
	auto checkPlanePlaneIntersection =
		[this](std::vector<MathLib::Point3d> & pnts) -> bool
	{
		// given intersection line
		MathLib::Vector3 d0(pnts[0]);
		// remove z component to obtain a horizontal line
		d0[2] = 0.0;
		MathLib::Vector3 const p0(pnts[1]);

		// construct arbitrary horizontal plane in Hessian normal form going through p0
		MathLib::Vector3 second_vector(pnts[2]);
		// again, as in d0, set z component to zero to obtain a horizontal plane
		second_vector[2] = 0.0;

		// Both planes go through p0; First plane is vertical. Second plane's
		// second spanning vector lies in horizontal plane.
		return this->check(
		    {d0, {0.0, 0.0, 1.0}, p0}, {d0, second_vector, p0}, d0);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(3,list_of(points_gen))),
		gtest_reporter);
}

TEST_F(GeoLibComputePlanePlaneIntersection,
	TestHorizontalPlaneXZPlaneIntersection)
{
	auto checkPlanePlaneIntersection =
		[this](std::vector<MathLib::Point3d> & pnts) -> bool
	{
		// create intersection line as a starting point for the test
		MathLib::Vector3 d0(pnts[0]);
		// remove y and z components -> line in parallel to the x axis
		d0[1] = 0.0;
		d0[2] = 0.0;
		MathLib::Vector3 const p0(pnts[1]);

		// construct arbitrary horizontal plane in Hessian normal form going through p0
		MathLib::Vector3 second_vector(pnts[2]);
		// set z component to zero to obtain a horizontal plane
		second_vector[2] = 0.0;

		// Both planes go through p0; First plane is vertical. Second plane's
		// second spanning vector lies in horizontal plane.
		return this->check(
		    {d0, {0.0, 0.0, 1.0}, p0}, {d0, second_vector, p0}, d0);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(3,list_of(points_gen))),
		gtest_reporter);
}

TEST_F(GeoLibComputePlanePlaneIntersection,
	TestHorizontalPlaneYZPlaneIntersection)
{
	auto checkPlanePlaneIntersection =
		[this](std::vector<MathLib::Point3d> & pnts) -> bool
	{
		// create intersection line as a starting point for the test
		MathLib::Vector3 d0(pnts[0]);
		// remove x and z components -> line in parallel to the y axis
		d0[0] = 0.0;
		d0[2] = 0.0;
		MathLib::Vector3 const p0(pnts[1]);

		// construct horizontal plane in Hessian normal form going through p0
		MathLib::Vector3 second_vector(pnts[2]);
		// set z component to zero to obtain a horizontal plane
		second_vector[2] = 0.0;

		// Both planes go through p0; First plane is vertical. Second plane's
		// second spanning vector lies in horizontal plane.
		return this->check(
		    {d0, {0.0, 0.0, 1.0}, p0}, {d0, second_vector, p0}, d0);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(3,list_of(points_gen))),
		gtest_reporter);
}

