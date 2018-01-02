/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>
#include <numeric>
#include "Tests/GeoLib/AutoCheckGenerators.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/LineSegment.h"

namespace ac = autocheck;

class GeoLibSortLineSegments : public testing::Test
{
public:
    using PointGenerator = ac::RandomCirclePointGeneratorXY<ac::generator<double>>;
    using SymmSegmentGenerator = ac::SymmSegmentGeneratorXY<PointGenerator>;

    PointGenerator point_generator = PointGenerator(
        MathLib::Point3d(std::array<double, 3>{{0.0, 0.0, 0.0}}), 1.0);
    SymmSegmentGenerator segment_generator = SymmSegmentGenerator{point_generator,
        std::bind(ac::reflect, point_generator.center, std::placeholders::_1)};

    ac::gtest_reporter gtest_reporter;
};

#if !defined(_MSC_VER) || (_MSC_VER >= 2000)
// Compilers of MVS below 2015 do not support unrestricted unions. The
// unrestricted union is used by autocheck to handle test data. The autocheck
// workaround for MVS compilers (below version 2015) contains a bug and in the
// consequence the tests crashes. For this reason the tests are disabled under
// this environments.

// Use a chord of the unit circle as the original line segment. The line segment
// will be partitioned into several sub segments. The set of subsegments are
// given to the algorithm.
TEST_F(GeoLibSortLineSegments, SortSubSegments)
{
    auto partitionSegment = [](GeoLib::LineSegment const& s0,
                               std::vector<std::size_t> const& sub_seg_ids,
                               double const dt,
                               std::vector<GeoLib::LineSegment>& sub_segments)
    {
        for (auto sub_seg_id : sub_seg_ids)
        {
            double t(dt * sub_seg_id);
            auto* sub_seg_begin_pnt(new GeoLib::Point{
                (1 - t) * s0.getBeginPoint()[0] + t * s0.getEndPoint()[0],
                (1 - t) * s0.getBeginPoint()[1] + t * s0.getEndPoint()[1],
                (1 - t) * s0.getBeginPoint()[2] + t * s0.getEndPoint()[2]});
            t += dt;
            auto* sub_seg_end_pnt(new GeoLib::Point{
                (1 - t) * s0.getBeginPoint()[0] + t * s0.getEndPoint()[0],
                (1 - t) * s0.getBeginPoint()[1] + t * s0.getEndPoint()[1],
                (1 - t) * s0.getBeginPoint()[2] + t * s0.getEndPoint()[2]});
            sub_segments.emplace_back(
                GeoLib::LineSegment{sub_seg_begin_pnt, sub_seg_end_pnt, true});
        }
    };

    auto checkSortedSubSegments = [](
        GeoLib::LineSegment const& s0,
        std::vector<GeoLib::LineSegment> const& sub_segments)
    {
        double eps(std::numeric_limits<double>::epsilon());
        if (MathLib::sqrDist(s0.getBeginPoint(),
                             sub_segments.front().getBeginPoint()) >= eps)
            return false;
        if (MathLib::sqrDist(s0.getEndPoint(),
                             sub_segments.back().getEndPoint()) >= eps)
            return false;
        for (std::size_t k(0); k < sub_segments.size() - 1; ++k)
        {
            if (MathLib::sqrDist(sub_segments[k].getEndPoint(),
                                 sub_segments[k + 1].getBeginPoint()) >= eps)
                return false;
        }
        return true;
    };

    auto testSortSegments = [partitionSegment, checkSortedSubSegments](
        GeoLib::LineSegment const& s0)
    {
        std::size_t const n_sub_segments(4);
        double const dt(1.0 / static_cast<double>(n_sub_segments));
        std::vector<std::size_t> sub_seg_ids(n_sub_segments);
        std::iota(sub_seg_ids.begin(), sub_seg_ids.end(), 0);
        do
        {
            std::vector<GeoLib::LineSegment> sub_segments;
            partitionSegment(s0, sub_seg_ids, dt, sub_segments);
            GeoLib::sortSegments(s0.getBeginPoint(), sub_segments);
            if (!checkSortedSubSegments(s0, sub_segments)) return false;
        } while (std::next_permutation(sub_seg_ids.begin(), sub_seg_ids.end()));
        return true;
    };

    ac::check<GeoLib::LineSegment>(
        testSortSegments, 100,
        ac::make_arbitrary(segment_generator),
        gtest_reporter);
}

#endif
