/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */


#include "LineSegment.h"

namespace GeoLib
{
LineSegment::LineSegment(Point* const a, Point* const b,
                         bool point_mem_management_by_line_segment)
    : _a(a),
      _b(b),
      _point_mem_management_by_line_segment(
          point_mem_management_by_line_segment)
{}

LineSegment::LineSegment(LineSegment const& line_segment)
    : _a(new Point(line_segment.getBeginPoint())),
      _b(new Point(line_segment.getEndPoint())),
      _point_mem_management_by_line_segment(true)
{}

LineSegment::LineSegment(LineSegment&& line_segment)
    : _a(line_segment._a),
      _b(line_segment._b),
      _point_mem_management_by_line_segment(
          line_segment._point_mem_management_by_line_segment)
{
    line_segment._a = nullptr;
    line_segment._b = nullptr;
    line_segment._point_mem_management_by_line_segment = false;
}

LineSegment::~LineSegment()
{
    if (_point_mem_management_by_line_segment) {
        delete _b;
        delete _a;
    }
}

LineSegment& LineSegment::operator=(LineSegment const&) = default;

LineSegment& LineSegment::operator=(LineSegment&& line_segment)
{
    _a = line_segment._a;
    _b = line_segment._b;
    _point_mem_management_by_line_segment =
        line_segment._point_mem_management_by_line_segment;

    line_segment._a = nullptr;
    line_segment._b = nullptr;
    line_segment._point_mem_management_by_line_segment = false;

    return *this;
}

Point const& LineSegment::getBeginPoint() const
{
    return *_a;
}

Point & LineSegment::getBeginPoint()
{
    return *_a;
}

Point const& LineSegment::getEndPoint() const
{
    return *_b;
}

Point & LineSegment::getEndPoint()
{
    return *_b;
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
}
