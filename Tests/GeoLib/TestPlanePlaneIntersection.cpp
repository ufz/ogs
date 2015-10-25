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

#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>

#include "Tests/MathLib/AutoCheckTools.h"

#include "MathLib/MathTools.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/Plane.h"

namespace ac = autocheck;

// Parametrize over the generator type for the d0, u, and v vectors.
template <typename GeneratorTypes>
struct ConstructedLine : public ::testing::Test
{
	// Select vector generator from the GeneratorTypes vector.
	ac::cons_generator<
	    MathLib::Vector3,
	    typename boost::mpl::at<GeneratorTypes, boost::mpl::int_<0>>::type>
	    d0_generator;
	ac::cons_generator<
	    MathLib::Vector3,
	    typename boost::mpl::at<GeneratorTypes, boost::mpl::int_<1>>::type>
	    u_generator;
	ac::cons_generator<
	    MathLib::Vector3,
	    typename boost::mpl::at<GeneratorTypes, boost::mpl::int_<2>>::type>
	    v_generator;

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
		    MathLib::Point3d const&,     // Common plane point
		    MathLib::Vector3 const& u,   // First plane's second spanning vector
		    MathLib::Vector3 const& v)  // Second plane's second spanning vector
		            {
			            // triple with entries, " ", "xy", "yz", "zx"
			            return "[" + classifyVector(d0) + ", " +
			                   classifyVector(u) + ", " + classifyVector(v) +
			                   "]";

			        });
	}
	ac::classifier<MathLib::Vector3,
	          MathLib::Point3d,
	          MathLib::Vector3,
	          MathLib::Vector3> cls;

	// Check correctness of the computePlanePlaneIntersection algorithm by
	// construction of two planes intersecting a line through point p0 in
	// direction d0 and comparing the result of the algorithm. The direction
	// vector must be parallel and the point p on the computed intersection line
	// must lie in both planes.
	static bool check(
	    MathLib::Vector3 const& d0,  // First spanning vector
	    MathLib::Point3d const& p0,  // Common plane point
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

TYPED_TEST_CASE_P(ConstructedLine);

TYPED_TEST_P(ConstructedLine, TestPlanePlaneIntersection)
{
	ac::check<MathLib::Vector3, MathLib::Point3d, MathLib::Vector3,
	          MathLib::Vector3>(
	    this->check, 1000,
	    ac::make_arbitrary(this->d0_generator, this->points_gen,
	                       this->u_generator, this->v_generator)
	        .discard_if([](MathLib::Vector3 const& d0,
	                       MathLib::Point3d const&,
	                       MathLib::Vector3 const& u,
	                       MathLib::Vector3 const& v)
	                    {
		                    MathLib::Vector3 const zero{0, 0, 0};
		                    return (d0 == zero || u == zero || v == zero);
		                }),
	    this->gtest_reporter, this->cls);
}

REGISTER_TYPED_TEST_CASE_P(ConstructedLine, TestPlanePlaneIntersection);

using R = ac::randomTupleGenerator<double, 3>;
using XY = ac::tripleInPlaneGenerator<ac::CartesianPlane::XY, double>;
using YZ = ac::tripleInPlaneGenerator<ac::CartesianPlane::YZ, double>;
using ZX = ac::tripleInPlaneGenerator<ac::CartesianPlane::ZX, double>;
using X = ac::tripleOnAxisGenerator<ac::CartesianAxes::X, double>;
using Y = ac::tripleOnAxisGenerator<ac::CartesianAxes::Y, double>;
using Z = ac::tripleOnAxisGenerator<ac::CartesianAxes::Z, double>;

// There are 7*7*7 possible cases, let's try few of the combinations.
typedef ::testing::Types<
    boost::mpl::vector<R ,R ,R >,
    boost::mpl::vector<R ,Z ,R >,
    boost::mpl::vector<XY,Z ,XY>,
    boost::mpl::vector<X ,Z ,XY>,
    boost::mpl::vector<Y ,Z ,XY>,

    boost::mpl::vector<X ,R ,R >,
    boost::mpl::vector<Y ,R ,R >,
    boost::mpl::vector<Z ,R ,R >,
    boost::mpl::vector<XY,R ,R >,
    boost::mpl::vector<YZ,R ,R >,
    boost::mpl::vector<ZX,R ,R >,

    boost::mpl::vector<R ,X ,R >,
    boost::mpl::vector<R ,Y ,R >,
    boost::mpl::vector<R ,Z ,R >,
    boost::mpl::vector<R ,XY,R >,
    boost::mpl::vector<R ,YZ,R >,
    boost::mpl::vector<R ,ZX,R >,

    boost::mpl::vector<R ,R ,X  >,
    boost::mpl::vector<R ,R ,Y  >,
    boost::mpl::vector<R ,R ,Z  >,
    boost::mpl::vector<R ,R ,XY >,
    boost::mpl::vector<R ,R ,YZ >,
    boost::mpl::vector<R ,R ,ZX >
> GeneratorTypes;

INSTANTIATE_TYPED_TEST_CASE_P(GeoLibComputePlanePlaneIntersection,
                              ConstructedLine, GeneratorTypes);
