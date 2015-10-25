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

	constexpr static double eps = std::numeric_limits<double>::epsilon();
	bool inXY(MathLib::Vector3 const& v) const { return std::abs(v[2]) < eps;}
	bool inYZ(MathLib::Vector3 const& v) const { return std::abs(v[0]) < eps; }
	bool inZX(MathLib::Vector3 const& v) const { return std::abs(v[1]) < eps; }
	bool inX(MathLib::Vector3 const& v) const { return inXY(v) && inZX(v); }
	bool inY(MathLib::Vector3 const& v) const { return inXY(v) && inYZ(v); }
	bool inZ(MathLib::Vector3 const& v) const { return inYZ(v) && inZX(v); }
	bool in0(MathLib::Vector3 const& v) const { return inX(v) && inYZ(v); }

	std::string classifyVector(MathLib::Vector3 const& v) const
	{
		if (in0(v)) return "0 ";
		if (inX(v)) return "x ";
		if (inY(v)) return "y ";
		if (inZ(v)) return "z ";
		if (inXY(v)) return "xy";
		if (inYZ(v)) return "yz";
		if (inZX(v)) return "zx";
		return "  ";	// not on any of the cartesian axes or planes.
	}

	void SetUp()
	{
		cls.collect([this](
		    MathLib::Vector3 const& d0,  // First spanning vector
		    MathLib::Vector3 const&,     // Common plane point
		    MathLib::Vector3 const& u,   // First plane's second spanning vector
		    MathLib::Vector3 const& v)  // Second plane's second spanning vector
		            {
			            // triple with entries, " ", "xy", "yz", "zx"
			            return "[" + classifyVector(d0) + ", " +
			                   classifyVector(u) + ", " + classifyVector(v) +
			                   "]";

			        });
	}
	ac::classifier<MathLib::Point3d,
	          MathLib::Point3d,
	          MathLib::Point3d,
	          MathLib::Point3d> cls;

	// Check correctness of the computePlanePlaneIntersection algorithm by
	// construction of two planes intersecting a line through point p0 in
	// direction d0 and comparing the result of the algorithm. The direction
	// vector must be parallel and the point p on the computed intersection line
	// must lie in both planes.
	static bool check(
	    MathLib::Vector3 const& d0,  // First spanning vector
	    MathLib::Vector3 const& p0,  // Common plane point
	    MathLib::Vector3 const& u,   // First plane's second spanning vector
	    MathLib::Vector3 const& v)   // Second plane's second spanning vector
	{
		GeoLib::Plane const plane1{d0, u, p0};
		GeoLib::Plane const plane2{d0, v, p0};
		// reconstructed intersection line
		MathLib::Vector3 d;
		MathLib::Point3d p;

		std::tie(d, p) = GeoLib::computePlanePlaneIntersection(plane1, plane2);

		// check if the given vector d0 and the computed vector are parallel
		if (!GeoLib::isParallel(d, d0))
		{
			std::cerr << "Expected d and d0 to be parallel.\n";
			return false;
		}

		// the point of the intersection line must lie inside both planes
		if (!plane1.isPointInPlane(p))
		{
			std::cerr << "Expected p to lie in plane1.\n";
			return false;
		}
		if (!plane2.isPointInPlane(p))
		{
			std::cerr << "Expected p to lie in plane1.\n";
			return false;
		}
		return true;
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
		return this->check(d0, p0, pnts[2], pnts[3]);
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
		return this->check(d0, p0, {0.0, 0.0, 1.0}, pnts[2]);
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
		return this->check(d0, p0, {0.0, 0.0, 1.0}, second_vector);
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
		return this->check(d0, p0, {0.0, 0.0, 1.0}, second_vector);
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
		return this->check(d0, p0, {0.0, 0.0, 1.0}, second_vector);
	};

	ac::check<std::vector<MathLib::Point3d>>(
		checkPlanePlaneIntersection,
		1000,
		ac::make_arbitrary(ac::fix(3,list_of(points_gen))),
		gtest_reporter);
}

