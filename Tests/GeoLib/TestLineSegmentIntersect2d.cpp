/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"
#include <array>
#include <ctime>
#include <functional>
#include <memory>
#include <random>

#include "Tests/GeoLib/AutoCheckGenerators.h"

#include "MathLib/Point3d.h"
#include "GeoLib/AnalyticalGeometry.h"

namespace ac = autocheck;

class LineSegmentIntersect2dTest : public testing::Test
{
public:
    using PointGenerator = ac::RandomCirclePointGeneratorXY<ac::generator<double>>;
    using SymmSegmentGenerator = ac::SymmSegmentGeneratorXY<PointGenerator>;
    using PairSegmentGenerator = ac::PairSegmentGeneratorXY<SymmSegmentGenerator>;

    PointGenerator point_generator1 = PointGenerator(
        MathLib::Point3d(std::array<double, 3>{{0.0, 0.0, 0.0}}), 1.0);
    SymmSegmentGenerator segment_generator1 = SymmSegmentGenerator{point_generator1,
        std::bind(ac::reflect, point_generator1.center, std::placeholders::_1)};

    PointGenerator point_generator2 = PointGenerator(
        MathLib::Point3d(std::array<double, 3>{{2.0, 0.0, 0.0}}), 1.0);
    SymmSegmentGenerator segment_generator2 = SymmSegmentGenerator{point_generator2,
        std::bind(ac::reflect, point_generator2.center, std::placeholders::_1)};

    MathLib::Vector3 const translation_vector1 = {2,2,0};
    PairSegmentGenerator pair_segment_generator1 = PairSegmentGenerator{
        segment_generator1,
        std::bind(ac::translate, translation_vector1, std::placeholders::_1)};

    MathLib::Vector3 const translation_vector2 = {0,0,0};
    PairSegmentGenerator pair_segment_generator2 = PairSegmentGenerator{
        segment_generator1,
        std::bind(ac::translate, translation_vector2, std::placeholders::_1)};

    ac::gtest_reporter gtest_reporter;
};

#if !defined(_MSC_VER) || (_MSC_VER >= 2000)
// Compilers of MVS below 2015 do not support unrestricted unions. The
// unrestricted union is used by autocheck to handle test data. The autocheck
// workaround for MVS compilers (below version 2015) contains a bug and in the
// consequence the tests crashes. For this reason the tests are disabled under
// this environments.

// Test the intersection of intersecting line segments. Line segments are chords
// of the same circle that both contains the center of the circle. As a
// consequence the center of the circle is the intersection point.
TEST_F(LineSegmentIntersect2dTest, RandomSegmentOrientationIntersecting)
{
    auto intersect = [](GeoLib::LineSegment const& s0,
                        GeoLib::LineSegment const& s1)
    {
        auto ipnts = GeoLib::lineSegmentIntersect2d(s0, s1);
        if (ipnts.size() == 1)
        {
            MathLib::Point3d const center{std::array<double, 3>{
                {(s0.getBeginPoint()[0] + s0.getEndPoint()[0]) / 2,
                 (s0.getBeginPoint()[1] + s0.getEndPoint()[1]) / 2, 0.0}}};
            const double sqr_dist(MathLib::sqrDist(ipnts[0], center));
            if (sqr_dist < std::numeric_limits<double>::epsilon())
                return true;
            else
                return false;
        }
        return ipnts.size() == 2;
    };

    ac::check<GeoLib::LineSegment, GeoLib::LineSegment>(
        intersect, 1000,
        ac::make_arbitrary(segment_generator1, segment_generator1),
        gtest_reporter);
}

// Test the intersection of non-intersecting line segments. Line segments are
// chords of non-intersecting circles.
TEST_F(LineSegmentIntersect2dTest, RandomSegmentOrientationNonIntersecting)
{
    auto intersect = [](GeoLib::LineSegment const& s0,
                        GeoLib::LineSegment const& s1)
    {
        auto ipnts = GeoLib::lineSegmentIntersect2d(s0, s1);
        return ipnts.empty();
    };

    // generate non-intersecting segments
    ac::check<GeoLib::LineSegment, GeoLib::LineSegment>(
        intersect, 1000,
        ac::make_arbitrary(segment_generator1, segment_generator2),
        gtest_reporter);
}

// Test the intersection of non-intersecting, parallel line segments. The second
// line segment is created by translating the first line segment.
TEST_F(LineSegmentIntersect2dTest, ParallelNonIntersectingSegmentOrientation)
{
    auto intersect = [](
        std::pair<GeoLib::LineSegment const&, GeoLib::LineSegment const&> const&
            segment_pair)
    {
        auto ipnts = GeoLib::lineSegmentIntersect2d(segment_pair.first,
                                                    segment_pair.second);
        return ipnts.empty();
    };

    // generate non-intersecting segments
    ac::check<std::pair<GeoLib::LineSegment, GeoLib::LineSegment>>(
        intersect, 1000,
        ac::make_arbitrary(pair_segment_generator1),
        gtest_reporter);
}

// Test the intersection of parallel, interfering line segments.
TEST_F(LineSegmentIntersect2dTest, ParallelIntersectingSegmentOrientation)
{
    auto intersect = [](
        std::pair<GeoLib::LineSegment const&, GeoLib::LineSegment const&> const&
            segment_pair)
    {
        auto ipnts = GeoLib::lineSegmentIntersect2d(segment_pair.first,
                                                    segment_pair.second);
        return ipnts.size() == 2;
    };

    // generate non-intersecting segments
    ac::check<std::pair<GeoLib::LineSegment, GeoLib::LineSegment>>(
        intersect, 1000,
        ac::make_arbitrary(pair_segment_generator2),
        gtest_reporter);
}

#endif
