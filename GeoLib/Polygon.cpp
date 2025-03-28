/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Implementation of the Polygon class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Polygon.h"

#include <algorithm>

#include "AnalyticalGeometry.h"
#include "BaseLib/quicksort.h"

namespace GeoLib
{
enum class Location
{
    LEFT,
    RIGHT,
    BEYOND,
    BEHIND,
    BETWEEN,
    SOURCE,
    DESTINATION
};

/**
 * edge classification
 */
enum class EdgeType
{
    TOUCHING,    //!< TOUCHING
    CROSSING,    //!< CROSSING
    INESSENTIAL  //!< INESSENTIAL
};

/**
 * 2D method - ignores z coordinate. It calculates the location
 * of the point relative to the line segment specified by the points source and
 * destination. (literature reference:
 * Computational Geometry and Computer Graphics in C++; Michael J. Laszlo)
 * \param source the first point of the line segment
 * \param destination the end point of the line segment
 * \param pnt the test point
 * \return a value of enum Location
 */
Location getLocationOfPoint(MathLib::Point3d const& source,
                            MathLib::Point3d const& destination,
                            MathLib::Point3d const& pnt)
{
    long double const a[2] = {destination[0] - source[0],
                              destination[1] - source[1]};  // vector
    long double const b[2] = {pnt[0] - source[0],
                              pnt[1] - source[1]};  // vector

    long double const det_2x2(a[0] * b[1] - a[1] * b[0]);
    constexpr double eps = std::numeric_limits<double>::epsilon();

    if (det_2x2 > eps)
    {
        return Location::LEFT;
    }
    if (eps < std::abs(det_2x2))
    {
        return Location::RIGHT;
    }
    if (a[0] * b[0] < 0.0 || a[1] * b[1] < 0.0)
    {
        return Location::BEHIND;
    }
    if (a[0] * a[0] + a[1] * a[1] < b[0] * b[0] + b[1] * b[1])
    {
        return Location::BEYOND;
    }
    if (MathLib::sqrDist(pnt, source) < pow(eps, 2))
    {
        return Location::SOURCE;
    }
    if (MathLib::sqrDist(pnt, destination) < std::sqrt(eps))
    {
        return Location::DESTINATION;
    }
    return Location::BETWEEN;
}

/**
 * from book: Computational Geometry and Computer Graphics in C++, page 119
 * get the type of edge with respect to the given point (2d method!)
 * \param a first point of line segment
 * \param b last point of line segment
 * \param pnt point that is edge type computed for
 * \return a value of enum EdgeType
 */
EdgeType getEdgeType(MathLib::Point3d const& a,
                     MathLib::Point3d const& b,
                     MathLib::Point3d const& pnt)
{
    switch (getLocationOfPoint(a, b, pnt))
    {
        case Location::LEFT:
        {
            if (a[1] < pnt[1] && pnt[1] <= b[1])
            {
                return EdgeType::CROSSING;
            }

            return EdgeType::INESSENTIAL;
        }
        case Location::RIGHT:
        {
            if (b[1] < pnt[1] && pnt[1] <= a[1])
            {
                return EdgeType::CROSSING;
            }

            return EdgeType::INESSENTIAL;
        }
        case Location::BETWEEN:
        case Location::SOURCE:
        case Location::DESTINATION:
            return EdgeType::TOUCHING;
        default:
            return EdgeType::INESSENTIAL;
    }
}

Polygon::Polygon(const Polyline& ply, bool init)
    : Polyline(ply), _aabb(ply.getPointsVec(), ply.getPolylinePointIDs())
{
    if (init)
    {
        initialise();
    }
    _simple_polygon_list.push_back(this);
}

Polygon::Polygon(Polygon const& other) : Polyline(other), _aabb(other._aabb)
{
    _simple_polygon_list.push_back(this);
    auto sub_polygon_it(other._simple_polygon_list.begin());
    for (sub_polygon_it++;  // the first entry is the polygon itself, skip the
                            // entry
         sub_polygon_it != other._simple_polygon_list.end();
         ++sub_polygon_it)
    {
        _simple_polygon_list.emplace_back(new Polygon(*(*sub_polygon_it)));
    }
}

Polygon::~Polygon()
{
    // remove polygons from list
    for (auto& polygon : _simple_polygon_list)
    {
        // the first entry of the list can be a pointer the object itself
        if (polygon != this)
        {
            delete polygon;
        }
    }
}

bool Polygon::initialise()
{
    if (this->isClosed())
    {
        ensureCCWOrientation();
        return true;
    }
    WARN("Polygon::initialise(): base polyline is not closed.");
    return false;
}

/**
 * Computes all intersections of the straight line segment and the polyline
 * boundary
 * \param polygon the polygon the segment line segment that will be processed
 * \param segment the line segment that will be processed
 * \return a possible empty vector containing the intersection points
 */
std::vector<GeoLib::Point> getAllIntersectionPoints(
    Polygon const& polygon, GeoLib::LineSegment const& segment)
{
    std::vector<GeoLib::Point> intersections;
    GeoLib::Point s;
    for (auto&& seg_it : polygon)
    {
        if (GeoLib::lineSegmentIntersect(seg_it, segment, s))
        {
            intersections.push_back(s);
        }
    }

    return intersections;
}

bool Polygon::isPntInPolygon(MathLib::Point3d const& pnt) const
{
    auto const [min_aabb_pnt, max_aabb_pnt] = _aabb.getMinMaxPoints();

    if (pnt[0] < min_aabb_pnt[0] || max_aabb_pnt[0] < pnt[0] ||
        pnt[1] < min_aabb_pnt[1] || max_aabb_pnt[1] < pnt[1])
    {
        return false;
    }

    if (_simple_polygon_list.size() == 1)
    {
        std::size_t n_intersections(0);
        const std::size_t n_nodes(getNumberOfPoints() - 1);
        for (std::size_t k(0); k < n_nodes; k++)
        {
            if (((*(getPoint(k)))[1] <= pnt[1] &&
                 pnt[1] <= (*(getPoint(k + 1)))[1]) ||
                ((*(getPoint(k + 1)))[1] <= pnt[1] &&
                 pnt[1] <= (*(getPoint(k)))[1]))
            {
                switch (getEdgeType(*getPoint(k), *getPoint(k + 1), pnt))
                {
                    case EdgeType::TOUCHING:
                        return true;
                    case EdgeType::CROSSING:
                        n_intersections++;
                        break;
                    case EdgeType::INESSENTIAL:
                        break;
                    default:
                        // do nothing
                        ;
                }
            }
        }
        if (n_intersections % 2 == 1)
        {
            return true;
        }
    }
    else
    {
        for (auto it(_simple_polygon_list.begin()++);
             it != _simple_polygon_list.end();
             ++it)
        {
            if ((*it)->isPntInPolygon(pnt))
            {
                return true;
            }
        }
    }
    return false;
}

bool Polygon::containsSegment(GeoLib::LineSegment const& segment) const
{
    std::vector<GeoLib::Point> s(getAllIntersectionPoints(*this, segment));

    GeoLib::Point const& a{segment.getBeginPoint()};
    GeoLib::Point const& b{segment.getEndPoint()};
    // no intersections -> check if at least one point of segment is in polygon
    if (s.empty())
    {
        return (isPntInPolygon(a));
    }

    const double tol(std::numeric_limits<float>::epsilon());

    // one intersection, intersection in line segment end point
    if (s.size() == 1)
    {
        const double sqr_dist_as(MathLib::sqrDist(a, s[0]));
        if (sqr_dist_as < tol)
        {
            return (isPntInPolygon(b));
        }

        const double sqr_dist_bs(MathLib::sqrDist(b, s[0]));
        if (sqr_dist_bs < tol)
        {
            return (isPntInPolygon(a));
        }
    }

    // Sorting the intersection with respect to the distance to the point a.
    // This induces a partition of the line segment into sub segments.
    std::sort(s.begin(), s.end(),
              [&a](GeoLib::Point const& p0, GeoLib::Point const& p1)
              { return MathLib::sqrDist(a, p0) < MathLib::sqrDist(a, p1); });

    // remove sub segments with almost zero length
    for (std::size_t k(0); k < s.size() - 1;)
    {
        if (MathLib::sqrDist(s[k], s[k + 1]) < tol)
        {
            s.erase(s.begin() + k + 1);
        }
        else
        {
            k++;
        }
    }

    // Check if all sub segments are within the polygon.
    if (!isPntInPolygon(GeoLib::Point(0.5 * (a[0] + s[0][0]),
                                      0.5 * (a[1] + s[0][1]),
                                      0.5 * (a[2] + s[0][2]))))
    {
        return false;
    }
    const std::size_t n_sub_segs(s.size() - 1);
    for (std::size_t k(0); k < n_sub_segs; k++)
    {
        if (!isPntInPolygon(GeoLib::Point(0.5 * (s[k][0] + s[k + 1][0]),
                                          0.5 * (s[k][1] + s[k + 1][1]),
                                          0.5 * (s[k][2] + s[k + 1][2]))))
        {
            return false;
        }
    }
    return isPntInPolygon(GeoLib::Point(0.5 * (s[0][0] + b[0]),
                                        0.5 * (s[0][1] + b[1]),
                                        0.5 * (s[0][2] + b[2])));
}

bool Polygon::isPolylineInPolygon(const Polyline& ply) const
{
    return std::all_of(ply.begin(), ply.end(),
                       [this](auto const& segment)
                       { return containsSegment(segment); });
}

bool Polygon::isPartOfPolylineInPolygon(const Polyline& ply) const
{
    const std::size_t ply_size(ply.getNumberOfPoints());
    // check points
    for (std::size_t k(0); k < ply_size; k++)
    {
        if (isPntInPolygon(*(ply.getPoint(k))))
        {
            return true;
        }
    }

    auto polygon_segment_intersects_line = [&](auto const& polygon_seg)
    {
        GeoLib::Point s;
        return std::any_of(ply.begin(), ply.end(),
                           [&polygon_seg, &s](auto const& polyline_seg) {
                               return GeoLib::lineSegmentIntersect(
                                   polyline_seg, polygon_seg, s);
                           });
    };

    return std::any_of(std::cbegin(*this), std::cend(*this),
                       polygon_segment_intersects_line);
}

bool Polygon::getNextIntersectionPointPolygonLine(
    GeoLib::LineSegment const& seg, GeoLib::Point& intersection_pnt,
    std::size_t& seg_num) const
{
    if (_simple_polygon_list.size() == 1)
    {
        for (auto seg_it(begin() + seg_num); seg_it != end(); ++seg_it)
        {
            if (GeoLib::lineSegmentIntersect(*seg_it, seg, intersection_pnt))
            {
                seg_num = seg_it.getSegmentNumber();
                return true;
            }
        }
    }
    else
    {
        for (auto const* polygon : _simple_polygon_list)
        {
            for (auto seg_it(polygon->begin()); seg_it != polygon->end();
                 ++seg_it)
            {
                if (GeoLib::lineSegmentIntersect(*seg_it, seg,
                                                 intersection_pnt))
                {
                    seg_num = seg_it.getSegmentNumber();
                    return true;
                }
            }
        }
    }
    return false;
}

void Polygon::ensureCCWOrientation()
{
    // *** pre processing: rotate points to xy-plan
    // *** copy points to vector - last point is identical to the first
    std::size_t n_pnts(this->getNumberOfPoints() - 1);
    std::vector<GeoLib::Point*> tmp_polygon_pnts;
    for (std::size_t k(0); k < n_pnts; k++)
    {
        tmp_polygon_pnts.push_back(new GeoLib::Point(*(this->getPoint(k))));
    }

    // rotate copied points into x-y-plane
    GeoLib::rotatePointsToXY(tmp_polygon_pnts);

    for (auto& tmp_polygon_pnt : tmp_polygon_pnts)
    {
        (*tmp_polygon_pnt)[2] =
            0.0;  // should be -= d but there are numerical errors
    }

    // *** get the left most upper point
    std::size_t min_x_max_y_idx(0);  // for orientation check
    for (std::size_t k(0); k < n_pnts; k++)
    {
        if ((*(tmp_polygon_pnts[k]))[0] <=
            (*(tmp_polygon_pnts[min_x_max_y_idx]))[0])
        {
            if ((*(tmp_polygon_pnts[k]))[0] <
                (*(tmp_polygon_pnts[min_x_max_y_idx]))[0])
            {
                min_x_max_y_idx = k;
            }
            else if ((*(tmp_polygon_pnts[k]))[1] >
                     (*(tmp_polygon_pnts[min_x_max_y_idx]))[1])
            {
                min_x_max_y_idx = k;
            }
        }
    }
    // *** determine orientation
    GeoLib::Orientation orient;
    if (0 < min_x_max_y_idx && min_x_max_y_idx < n_pnts - 2)
    {
        orient = GeoLib::getOrientation(*tmp_polygon_pnts[min_x_max_y_idx - 1],
                                        *tmp_polygon_pnts[min_x_max_y_idx],
                                        *tmp_polygon_pnts[min_x_max_y_idx + 1]);
    }
    else
    {
        if (0 == min_x_max_y_idx)
        {
            orient = GeoLib::getOrientation(*tmp_polygon_pnts[n_pnts - 1],
                                            *tmp_polygon_pnts[0],
                                            *tmp_polygon_pnts[1]);
        }
        else
        {
            orient = GeoLib::getOrientation(*tmp_polygon_pnts[n_pnts - 2],
                                            *tmp_polygon_pnts[n_pnts - 1],
                                            *tmp_polygon_pnts[0]);
        }
    }

    if (orient != GeoLib::CCW)
    {
        reverseOrientation();
    }

    for (std::size_t k(0); k < n_pnts; k++)
    {
        delete tmp_polygon_pnts[k];
    }
}

void Polygon::splitPolygonAtIntersection(
    const std::list<Polygon*>::const_iterator& polygon_it)
{
    GeoLib::Polyline::SegmentIterator seg_it0((*polygon_it)->begin());
    GeoLib::Polyline::SegmentIterator seg_it1((*polygon_it)->begin());
    GeoLib::Point intersection_pnt;
    if (!GeoLib::lineSegmentsIntersect(*polygon_it, seg_it0, seg_it1,
                                       intersection_pnt))
    {
        return;
    }

    std::size_t idx0(seg_it0.getSegmentNumber());
    std::size_t idx1(seg_it1.getSegmentNumber());
    // adding intersection point to pnt_vec
    std::size_t const intersection_pnt_id(_ply_pnts.size());
    const_cast<std::vector<Point*>&>(_ply_pnts).push_back(
        new GeoLib::Point(intersection_pnt));

    // split Polygon
    if (idx0 > idx1)
    {
        std::swap(idx0, idx1);
    }

    GeoLib::Polyline polyline0{(*polygon_it)->getPointsVec()};
    for (std::size_t k(0); k <= idx0; k++)
    {
        polyline0.addPoint((*polygon_it)->getPointID(k));
    }
    polyline0.addPoint(intersection_pnt_id);
    for (std::size_t k(idx1 + 1); k < (*polygon_it)->getNumberOfPoints(); k++)
    {
        polyline0.addPoint((*polygon_it)->getPointID(k));
    }

    GeoLib::Polyline polyline1{(*polygon_it)->getPointsVec()};
    polyline1.addPoint(intersection_pnt_id);
    for (std::size_t k(idx0 + 1); k <= idx1; k++)
    {
        polyline1.addPoint((*polygon_it)->getPointID(k));
    }
    polyline1.addPoint(intersection_pnt_id);

    // remove the polygon except the first
    if (*polygon_it != this)
    {
        delete *polygon_it;
    }
    // erase polygon_it and add two new polylines
    auto polygon1_it = _simple_polygon_list.insert(
        _simple_polygon_list.erase(polygon_it), new GeoLib::Polygon(polyline1));
    auto polygon0_it = _simple_polygon_list.insert(
        polygon1_it, new GeoLib::Polygon(polyline0));

    splitPolygonAtIntersection(polygon0_it);
    splitPolygonAtIntersection(polygon1_it);
}

void Polygon::splitPolygonAtPoint(
    const std::list<GeoLib::Polygon*>::iterator& polygon_it)
{
    std::size_t const n((*polygon_it)->getNumberOfPoints() - 1);
    std::vector<std::size_t> id_vec(n);
    std::vector<std::size_t> perm(n);
    for (std::size_t k(0); k < n; k++)
    {
        id_vec[k] = (*polygon_it)->getPointID(k);
        perm[k] = k;
    }

    BaseLib::quicksort(id_vec, 0, n, perm);

    for (std::size_t k(0); k < n - 1; k++)
    {
        if (id_vec[k] == id_vec[k + 1])
        {
            std::size_t idx0 = perm[k];
            std::size_t idx1 = perm[k + 1];

            if (idx0 > idx1)
            {
                std::swap(idx0, idx1);
            }

            // create two closed polylines
            GeoLib::Polyline polyline0{*(*polygon_it)};
            for (std::size_t j(0); j <= idx0; j++)
            {
                polyline0.addPoint((*polygon_it)->getPointID(j));
            }
            for (std::size_t j(idx1 + 1);
                 j < (*polygon_it)->getNumberOfPoints();
                 j++)
            {
                polyline0.addPoint((*polygon_it)->getPointID(j));
            }

            GeoLib::Polyline polyline1{*(*polygon_it)};
            for (std::size_t j(idx0); j <= idx1; j++)
            {
                polyline1.addPoint((*polygon_it)->getPointID(j));
            }

            // remove the polygon except the first
            if (*polygon_it != this)
            {
                delete *polygon_it;
            }
            // erase polygon_it and add two new polygons
            auto polygon1_it = _simple_polygon_list.insert(
                _simple_polygon_list.erase(polygon_it), new Polygon(polyline1));
            auto polygon0_it = _simple_polygon_list.insert(
                polygon1_it, new Polygon(polyline0));

            splitPolygonAtPoint(polygon0_it);
            splitPolygonAtPoint(polygon1_it);

            return;
        }
    }
}

bool operator==(Polygon const& lhs, Polygon const& rhs)
{
    if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
    {
        return false;
    }

    const std::size_t n(lhs.getNumberOfPoints());
    const std::size_t start_pnt(lhs.getPointID(0));

    // search start point of first polygon in second polygon
    bool nfound(true);
    std::size_t k(0);
    for (; k < n - 1 && nfound; k++)
    {
        if (start_pnt == rhs.getPointID(k))
        {
            nfound = false;
            break;
        }
    }

    // case: start point not found in second polygon
    if (nfound)
    {
        return false;
    }

    // *** determine direction
    // opposite direction
    if (k == n - 2)
    {
        for (k = 1; k < n - 1; k++)
        {
            if (lhs.getPointID(k) != rhs.getPointID(n - 1 - k))
            {
                return false;
            }
        }
        return true;
    }

    // same direction - start point of first polygon at arbitrary position in
    // second polygon
    if (lhs.getPointID(1) == rhs.getPointID(k + 1))
    {
        std::size_t j(k + 2);
        for (; j < n - 1; j++)
        {
            if (lhs.getPointID(j - k) != rhs.getPointID(j))
            {
                return false;
            }
        }
        j = 0;  // new start point at second polygon
        for (; j < k + 1; j++)
        {
            if (lhs.getPointID(n - (k + 2) + j + 1) != rhs.getPointID(j))
            {
                return false;
            }
        }
        return true;
    }
    // opposite direction with start point of first polygon at arbitrary
    // position
    // *** ATTENTION
    WARN(
        "operator==(Polygon const& lhs, Polygon const& rhs) - not tested case "
        "(implementation is probably buggy) - please contact "
        "thomas.fischer@ufz.de mentioning the problem.");
    // in second polygon
    if (lhs.getPointID(1) == rhs.getPointID(k - 1))
    {
        std::size_t j(k - 2);
        for (; j > 0; j--)
        {
            if (lhs.getPointID(k - 2 - j) != rhs.getPointID(j))
            {
                return false;
            }
        }
        // new start point at second polygon - the point n-1 of a polygon is
        // equal to the first point of the polygon (for this reason: n-2)
        j = n - 2;
        for (; j > k - 1; j--)
        {
            if (lhs.getPointID(n - 2 + j + k - 2) != rhs.getPointID(j))
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

std::list<Polygon*> const& Polygon::computeListOfSimplePolygons()
{
    splitPolygonAtPoint(_simple_polygon_list.begin());
    splitPolygonAtIntersection(_simple_polygon_list.begin());

    for (auto polygon : _simple_polygon_list)
    {
        polygon->initialise();
    }
    return _simple_polygon_list;
}

}  // end namespace GeoLib
