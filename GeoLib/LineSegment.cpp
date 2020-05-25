/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LineSegment.h"

namespace GeoLib
{
LineSegment::LineSegment(Point* const a, Point* const b,
                         bool point_mem_management_by_line_segment)
    : a_(a),
      b_(b),
      point_mem_management_by_line_segment_(
          point_mem_management_by_line_segment)
{}

LineSegment::LineSegment(LineSegment const& line_segment)
    : a_(new Point(line_segment.getBeginPoint())),
      b_(new Point(line_segment.getEndPoint())),
      point_mem_management_by_line_segment_(true)
{}

LineSegment::LineSegment(LineSegment&& line_segment)
    : a_(line_segment.a_),
      b_(line_segment.b_),
      point_mem_management_by_line_segment_(
          line_segment.point_mem_management_by_line_segment_)
{
    line_segment.a_ = nullptr;
    line_segment.b_ = nullptr;
    line_segment.point_mem_management_by_line_segment_ = false;
}

LineSegment::~LineSegment()
{
    if (point_mem_management_by_line_segment_) {
        delete b_;
        delete a_;
    }
}

LineSegment& LineSegment::operator=(LineSegment const&) = default;

LineSegment& LineSegment::operator=(LineSegment&& line_segment)
{
    a_ = line_segment.a_;
    b_ = line_segment.b_;
    point_mem_management_by_line_segment_ =
        line_segment.point_mem_management_by_line_segment_;

    line_segment.a_ = nullptr;
    line_segment.b_ = nullptr;
    line_segment.point_mem_management_by_line_segment_ = false;

    return *this;
}

Point const& LineSegment::getBeginPoint() const
{
    return *a_;
}

Point & LineSegment::getBeginPoint()
{
    return *a_;
}

Point const& LineSegment::getEndPoint() const
{
    return *b_;
}

Point & LineSegment::getEndPoint()
{
    return *b_;
}

std::ostream& operator<< (std::ostream& os, LineSegment const& s)
{
    os << "{(" << s.getBeginPoint() << "), (" << s.getEndPoint() << ")}";
    return os;
}

std::ostream& operator<<(std::ostream& os,
                         std::pair<GeoLib::LineSegment const&,
                                   GeoLib::LineSegment const&> const& seg_pair)
{
    os << seg_pair.first << " x " << seg_pair.second;
    return os;
}

bool operator==(LineSegment const& s0, LineSegment const& s1)
{
    double const tol(std::numeric_limits<double>::epsilon());
    return (MathLib::sqrDist(s0.getBeginPoint(), s1.getBeginPoint()) < tol &&
         MathLib::sqrDist(s0.getEndPoint(), s1.getEndPoint()) < tol) ||
        (MathLib::sqrDist(s0.getBeginPoint(), s1.getEndPoint()) < tol &&
         MathLib::sqrDist(s0.getEndPoint(), s1.getBeginPoint()) < tol);
}
}  // namespace GeoLib
