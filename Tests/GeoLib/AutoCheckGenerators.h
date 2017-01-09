/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <cmath>
#include <memory>
#include <random>

#include "autocheck/autocheck.hpp"
#include "MathLib/Point3d.h"
#include "MathLib/Vector3.h"
#include "GeoLib/LineSegment.h"

namespace autocheck
{

// generates points on a circle in x-y plane
template <typename Gen = generator<double>>
struct RandomCirclePointGeneratorXY
{
    RandomCirclePointGeneratorXY(
        MathLib::Point3d const& c = MathLib::Point3d{std::array<double, 3>{{0, 0, 0}}},
        double r = 1.0)
        : center(c), radius(r)
    {}

    using result_type = MathLib::Point3d;

    result_type operator()(std::size_t /*size*/ = 0)
    {
        // generates uniformly distributed values in the interval [-4, 4]
        auto const angle(generator(4));
        return MathLib::Point3d{
            std::array<double, 3>{{center[0] + radius * std::cos(angle),
                                   center[1] + radius * std::sin(angle), 0.0}}};
    }

    Gen generator;
    MathLib::Point3d const center;
    double const radius;
};

// reflect point p on the point c in x-y plane
MathLib::Point3d reflect(MathLib::Point3d const& c, MathLib::Point3d const& p);

// generates line segments that are special chords of a circle, i.e. the chords
// include the circle centre point.
template <typename Gen = RandomCirclePointGeneratorXY<generator<double>>>
struct SymmSegmentGeneratorXY
{
    SymmSegmentGeneratorXY(
        Gen s,
        std::function<MathLib::Point3d(MathLib::Point3d const&)> f)
        : source(s), function(f)
    {}

    using result_type = GeoLib::LineSegment;

    result_type operator()(std::size_t /*size*/ = 0)
    {
        std::unique_ptr<GeoLib::Point> a{new GeoLib::Point{source(), 0}};
        std::unique_ptr<GeoLib::Point> b{new GeoLib::Point{function(*a), 1}};
        return result_type{a.release(), b.release(), true};
    }

    Gen source;
    std::function<MathLib::Point3d(MathLib::Point3d const&)> function;
};

GeoLib::LineSegment translate(MathLib::Vector3 const& translation,
                              GeoLib::LineSegment const& line_seg);

template <typename Gen = SymmSegmentGeneratorXY<
              RandomCirclePointGeneratorXY<generator<double>>>>
struct PairSegmentGeneratorXY
{
    PairSegmentGeneratorXY(
        Gen s, std::function<GeoLib::LineSegment(GeoLib::LineSegment const&)> f)
        : segment_generator(s), function(f)
    {
    }

    using result_type = std::pair<GeoLib::LineSegment, GeoLib::LineSegment>;

    result_type operator()(std::size_t /*size*/ = 0)
    {
        auto first = segment_generator();
        auto second = function(first);
        return result_type{first, second};
    }

    Gen segment_generator;
    std::function<GeoLib::LineSegment(GeoLib::LineSegment const&)> function;
};

}  // namespace autocheck
