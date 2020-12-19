/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Implementation of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AnalyticalGeometry.h"

#include <algorithm>
#include <cmath>
#include <limits>


#include <Eigen/Dense>

#include "BaseLib/StringTools.h"

#include "Polyline.h"
#include "PointVec.h"

#include "MathLib/GeometricBasics.h"

extern double orient2d(double *, double *, double *);
extern double orient2dfast(double*, double*, double*);

namespace ExactPredicates
{
double getOrientation2d(MathLib::Point3d const& a,
    MathLib::Point3d const& b, MathLib::Point3d const& c)
{
    return orient2d(const_cast<double*>(a.getCoords()),
        const_cast<double*>(b.getCoords()),
        const_cast<double*>(c.getCoords()));
}

double getOrientation2dFast(MathLib::Point3d const& a,
                            MathLib::Point3d const& b,
                            MathLib::Point3d const& c)
{
    return orient2dfast(const_cast<double*>(a.getCoords()),
                        const_cast<double*>(b.getCoords()),
                        const_cast<double*>(c.getCoords()));
}
}  // namespace ExactPredicates

namespace GeoLib
{
Orientation getOrientation(MathLib::Point3d const& p0,
                           MathLib::Point3d const& p1,
                           MathLib::Point3d const& p2)
{
    double const orientation = ExactPredicates::getOrientation2d(p0, p1, p2);
    if (orientation > 0)
    {
        return CCW;
    }
    if (orientation < 0)
    {
        return CW;
    }
    return COLLINEAR;
}

Orientation getOrientationFast(MathLib::Point3d const& p0,
                               MathLib::Point3d const& p1,
                               MathLib::Point3d const& p2)
{
    double const orientation =
        ExactPredicates::getOrientation2dFast(p0, p1, p2);
    if (orientation > 0)
    {
        return CCW;
    }
    if (orientation < 0)
    {
        return CW;
    }
    return COLLINEAR;
}

bool parallel(Eigen::Vector3d v, Eigen::Vector3d w)
{
    const double eps(std::numeric_limits<double>::epsilon());
    double const eps_squared = eps * eps;

    // check degenerated cases
    if (v.squaredNorm() < eps_squared)
    {
        return false;
    }

    if (w.squaredNorm() < eps_squared)
    {
        return false;
    }

    v.normalize();
    w.normalize();

    bool parallel(true);
    if (std::abs(v[0] - w[0]) > eps)
    {
        parallel = false;
    }
    if (std::abs(v[1] - w[1]) > eps)
    {
        parallel = false;
    }
    if (std::abs(v[2] - w[2]) > eps)
    {
        parallel = false;
    }

    if (! parallel) {
        parallel = true;
        // change sense of direction of v_normalised
        v *= -1.0;
        // check again
        if (std::abs(v[0] - w[0]) > eps)
        {
            parallel = false;
        }
        if (std::abs(v[1] - w[1]) > eps)
        {
            parallel = false;
        }
        if (std::abs(v[2] - w[2]) > eps)
        {
            parallel = false;
        }
    }

    return parallel;
}

bool lineSegmentIntersect(GeoLib::LineSegment const& s0,
                          GeoLib::LineSegment const& s1,
                          GeoLib::Point& s)
{
    GeoLib::Point const& pa{s0.getBeginPoint()};
    GeoLib::Point const& pb{s0.getEndPoint()};
    GeoLib::Point const& pc{s1.getBeginPoint()};
    GeoLib::Point const& pd{s1.getEndPoint()};

    if (!isCoplanar(pa, pb, pc, pd))
    {
        return false;
    }

    auto const a =
        Eigen::Map<Eigen::Vector3d const>(s0.getBeginPoint().getCoords());
    auto const b =
        Eigen::Map<Eigen::Vector3d const>(s0.getEndPoint().getCoords());
    auto const c =
        Eigen::Map<Eigen::Vector3d const>(s1.getBeginPoint().getCoords());
    auto const d =
        Eigen::Map<Eigen::Vector3d const>(s1.getEndPoint().getCoords());

    Eigen::Vector3d const v = b - a;
    Eigen::Vector3d const w = d - c;
    Eigen::Vector3d const qp = c - a;
    Eigen::Vector3d const pq = a - c;

    double const eps = std::numeric_limits<double>::epsilon();
    double const squared_eps = eps * eps;
    // handle special cases here to avoid computing intersection numerical
    if (qp.squaredNorm() < squared_eps || (d - a).squaredNorm() < squared_eps)
    {
        s = pa;
        return true;
    }
    if ((c - b).squaredNorm() < squared_eps ||
        (d - b).squaredNorm() < squared_eps)
    {
        s = pb;
        return true;
    }

    auto isLineSegmentIntersectingAB = [&v](Eigen::Vector3d const& ap,
                                            std::size_t i)
    {
        // check if p is located at v=(a,b): (ap = t*v, t in [0,1])
        return 0.0 <= ap[i] / v[i] && ap[i] / v[i] <= 1.0;
    };

    if (parallel(v,w)) { // original line segments (a,b) and (c,d) are parallel
        if (parallel(pq,v)) { // line segment (a,b) and (a,c) are also parallel
            // Here it is already checked that the line segments (a,b) and (c,d)
            // are parallel. At this point it is also known that the line
            // segment (a,c) is also parallel to (a,b). In that case it is
            // possible to express c as c(t) = a + t * (b-a) (analog for the
            // point d). Since the evaluation of all three coordinate equations
            // (x,y,z) have to lead to the same solution for the parameter t it
            // is sufficient to evaluate t only once.

            // Search id of coordinate with largest absolute value which is will
            // be used in the subsequent computations. This prevents division by
            // zero in case the line segments are parallel to one of the
            // coordinate axis.
            std::size_t i_max(std::abs(v[0]) <= std::abs(v[1]) ? 1 : 0);
            i_max = std::abs(v[i_max]) <= std::abs(v[2]) ? 2 : i_max;
            if (isLineSegmentIntersectingAB(qp, i_max)) {
                s = pc;
                return true;
            }
            Eigen::Vector3d const ad = d - a;
            if (isLineSegmentIntersectingAB(ad, i_max))
            {
                s = pd;
                return true;
            }
            return false;
        }
        return false;
    }

    // general case
    const double sqr_len_v(v.squaredNorm());
    const double sqr_len_w(w.squaredNorm());

    Eigen::Matrix2d mat;
    mat(0,0) = sqr_len_v;
    mat(0,1) = -v.dot(w);
    mat(1,1) = sqr_len_w;
    mat(1,0) = mat(0,1);

    Eigen::Vector2d rhs{v.dot(qp), w.dot(pq)};

    rhs = mat.partialPivLu().solve(rhs);

    // no theory for the following tolerances, determined by testing
    // lower tolerance: little bit smaller than zero
    const double l(-1.0*std::numeric_limits<float>::epsilon());
    // upper tolerance a little bit greater than one
    const double u(1.0+std::numeric_limits<float>::epsilon());
    if (rhs[0] < l || u < rhs[0] || rhs[1] < l || u < rhs[1]) {
        return false;
    }

    // compute points along line segments with minimal distance
    GeoLib::Point const p0(a[0]+rhs[0]*v[0], a[1]+rhs[0]*v[1], a[2]+rhs[0]*v[2]);
    GeoLib::Point const p1(c[0]+rhs[1]*w[0], c[1]+rhs[1]*w[1], c[2]+rhs[1]*w[2]);

    double const min_dist(std::sqrt(MathLib::sqrDist(p0, p1)));
    double const min_seg_len(
        std::min(std::sqrt(sqr_len_v), std::sqrt(sqr_len_w)));
    if (min_dist < min_seg_len * 1e-6) {
        s[0] = 0.5 * (p0[0] + p1[0]);
        s[1] = 0.5 * (p0[1] + p1[1]);
        s[2] = 0.5 * (p0[2] + p1[2]);
        return true;
    }

    return false;
}

bool lineSegmentsIntersect(const GeoLib::Polyline* ply,
                           GeoLib::Polyline::SegmentIterator &seg_it0,
                           GeoLib::Polyline::SegmentIterator &seg_it1,
                           GeoLib::Point& intersection_pnt)
{
    std::size_t const n_segs(ply->getNumberOfSegments());
    // Neighbouring segments always intersects at a common vertex. The algorithm
    // checks for intersections of non-neighbouring segments.
    for (seg_it0 = ply->begin(); seg_it0 != ply->end() - 2; ++seg_it0)
    {
        seg_it1 = seg_it0+2;
        std::size_t const seg_num_0 = seg_it0.getSegmentNumber();
        for ( ; seg_it1 != ply->end(); ++seg_it1) {
            // Do not check first and last segment, because they are
            // neighboured.
            if (!(seg_num_0 == 0 && seg_it1.getSegmentNumber() == n_segs - 1)) {
                if (lineSegmentIntersect(*seg_it0, *seg_it1, intersection_pnt)) {
                    return true;
                }
            }
        }
    }
    return false;
}

void rotatePoints(Eigen::Matrix3d const& rot_mat,
                  std::vector<GeoLib::Point*>& pnts)
{
    rotatePoints(rot_mat, pnts.begin(), pnts.end());
}

Eigen::Matrix3d computeRotationMatrixToXY(Eigen::Vector3d const& n)
{
    Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Zero();
    // check if normal points already in the right direction
    if (n[0] == 0 && n[1] == 0)
    {
        rot_mat(1, 1) = 1.0;

        if (n[2] > 0)
        {
            // identity matrix
            rot_mat(0, 0) = 1.0;
            rot_mat(2, 2) = 1.0;
        }
        else
        {
            // rotate by pi about the y-axis
            rot_mat(0, 0) = -1.0;
            rot_mat(2, 2) = -1.0;
        }

        return rot_mat;
    }

    // sqrt (n_1^2 + n_3^2)
    double const h0(std::sqrt(n[0] * n[0] + n[2] * n[2]));

    // In case the x and z components of the normal are both zero the rotation
    // to the x-z-plane is not required, i.e. only the rotation in the z-axis is
    // required. The angle is either pi/2 or 3/2*pi. Thus the components of
    // rot_mat are as follows.
    if (h0 < std::numeric_limits<double>::epsilon())
    {
        rot_mat(0, 0) = 1.0;
        if (n[1] > 0)
        {
            rot_mat(1, 2) = -1.0;
            rot_mat(2, 1) = 1.0;
        }
        else
        {
            rot_mat(1, 2) = 1.0;
            rot_mat(2, 1) = -1.0;
        }
        return rot_mat;
    }

    double const h1(1 / n.norm());

    // general case: calculate entries of rotation matrix
    rot_mat(0, 0) = n[2] / h0;
    rot_mat(0, 1) = 0;
    rot_mat(0, 2) = -n[0] / h0;
    rot_mat(1, 0) = -n[1] * n[0] / h0 * h1;
    rot_mat(1, 1) = h0 * h1;
    rot_mat(1, 2) = -n[1] * n[2] / h0 * h1;
    rot_mat(2, 0) = n[0] * h1;
    rot_mat(2, 1) = n[1] * h1;
    rot_mat(2, 2) = n[2] * h1;

    return rot_mat;
}

Eigen::Matrix3d rotatePointsToXY(std::vector<GeoLib::Point*>& pnts)
{
    return rotatePointsToXY(pnts.begin(), pnts.end(), pnts.begin(), pnts.end());
}

std::unique_ptr<GeoLib::Point> triangleLineIntersection(
    MathLib::Point3d const& a, MathLib::Point3d const& b,
    MathLib::Point3d const& c, MathLib::Point3d const& p,
    MathLib::Point3d const& q)
{
    auto const va = Eigen::Map<Eigen::Vector3d const>(a.getCoords());
    auto const vb = Eigen::Map<Eigen::Vector3d const>(b.getCoords());
    auto const vc = Eigen::Map<Eigen::Vector3d const>(c.getCoords());
    auto const vp = Eigen::Map<Eigen::Vector3d const>(p.getCoords());
    auto const vq = Eigen::Map<Eigen::Vector3d const>(q.getCoords());

    Eigen::Vector3d const pq = vq - vp;
    Eigen::Vector3d const pa = va - vp;
    Eigen::Vector3d const pb = vb - vp;
    Eigen::Vector3d const pc = vc - vp;

    double u = pq.cross(pc).dot(pb);
    if (u < 0)
    {
        return nullptr;
    }
    double v = pq.cross(pa).dot(pc);
    if (v < 0)
    {
        return nullptr;
    }
    double w = pq.cross(pb).dot(pa);
    if (w < 0)
    {
        return nullptr;
    }

    const double denom(1.0 / (u + v + w));
    u *= denom;
    v *= denom;
    w *= denom;
    return std::make_unique<GeoLib::Point>(u * a[0] + v * b[0] + w * c[0],
                                           u * a[1] + v * b[1] + w * c[1],
                                           u * a[2] + v * b[2] + w * c[2]);
}

void computeAndInsertAllIntersectionPoints(GeoLib::PointVec &pnt_vec,
    std::vector<GeoLib::Polyline*> & plys)
{
    auto computeSegmentIntersections = [&pnt_vec](GeoLib::Polyline& poly0,
                                                  GeoLib::Polyline& poly1)
    {
        for (auto seg0_it(poly0.begin()); seg0_it != poly0.end(); ++seg0_it)
        {
            for (auto seg1_it(poly1.begin()); seg1_it != poly1.end(); ++seg1_it)
            {
                GeoLib::Point s(0.0, 0.0, 0.0, pnt_vec.size());
                if (lineSegmentIntersect(*seg0_it, *seg1_it, s))
                {
                    std::size_t const id(
                        pnt_vec.push_back(new GeoLib::Point(s)));
                    poly0.insertPoint(seg0_it.getSegmentNumber() + 1, id);
                    poly1.insertPoint(seg1_it.getSegmentNumber() + 1, id);
                }
            }
        }
    };

    for (auto it0(plys.begin()); it0 != plys.end(); ++it0) {
        auto it1(it0);
        ++it1;
        for (; it1 != plys.end(); ++it1) {
            computeSegmentIntersections(*(*it0), *(*it1));
        }
    }
}

GeoLib::Polygon rotatePolygonToXY(GeoLib::Polygon const& polygon_in,
                                  Eigen::Vector3d& plane_normal)
{
    // 1 copy all points
    auto* polygon_pnts(new std::vector<GeoLib::Point*>);
    for (std::size_t k(0); k < polygon_in.getNumberOfPoints(); k++)
    {
        polygon_pnts->push_back(new GeoLib::Point(*(polygon_in.getPoint(k))));
    }

    // 2 rotate points
    double d_polygon;
    std::tie(plane_normal, d_polygon) = GeoLib::getNewellPlane(*polygon_pnts);
    Eigen::Matrix3d const rot_mat =
        GeoLib::computeRotationMatrixToXY(plane_normal);
    GeoLib::rotatePoints(rot_mat, *polygon_pnts);

    // 3 set z coord to zero
    std::for_each(polygon_pnts->begin(), polygon_pnts->end(),
        [] (GeoLib::Point* p) { (*p)[2] = 0.0; }
    );

    // 4 create new polygon
    GeoLib::Polyline rot_polyline(*polygon_pnts);
    for (std::size_t k(0); k < polygon_in.getNumberOfPoints(); k++)
    {
        rot_polyline.addPoint(k);
    }
    rot_polyline.addPoint(0);
    return GeoLib::Polygon(rot_polyline);
}

std::vector<MathLib::Point3d> lineSegmentIntersect2d(
    GeoLib::LineSegment const& ab, GeoLib::LineSegment const& cd)
{
    GeoLib::Point const& a{ab.getBeginPoint()};
    GeoLib::Point const& b{ab.getEndPoint()};
    GeoLib::Point const& c{cd.getBeginPoint()};
    GeoLib::Point const& d{cd.getEndPoint()};

    double const orient_abc(getOrientation(a, b, c));
    double const orient_abd(getOrientation(a, b, d));

    // check if the segment (cd) lies on the left or on the right of (ab)
    if ((orient_abc > 0 && orient_abd > 0) || (orient_abc < 0 && orient_abd < 0)) {
        return std::vector<MathLib::Point3d>();
    }

    // check: (cd) and (ab) are on the same line
    if (orient_abc == 0.0 && orient_abd == 0.0) {
        double const eps(std::numeric_limits<double>::epsilon());
        if (MathLib::sqrDist2d(a, c) < eps && MathLib::sqrDist2d(b, d) < eps)
        {
            return {{a, b}};
        }
        if (MathLib::sqrDist2d(a, d) < eps && MathLib::sqrDist2d(b, c) < eps)
        {
            return {{a, b}};
        }

        // Since orient_ab and orient_abd vanish, a, b, c, d are on the same
        // line and for this reason it is enough to check the x-component.
        auto isPointOnSegment = [](double q, double p0, double p1)
        {
            double const t((q - p0) / (p1 - p0));
            return 0 <= t && t <= 1;
        };

        // check if c in (ab)
        if (isPointOnSegment(c[0], a[0], b[0])) {
            // check if a in (cd)
            if (isPointOnSegment(a[0], c[0], d[0])) {
                return {{a, c}};
            }
            // check b == c
            if (MathLib::sqrDist2d(b,c) < eps) {
                return {{b}};
            }
            // check if b in (cd)
            if (isPointOnSegment(b[0], c[0], d[0])) {
                return {{b, c}};
            }
            // check d in (ab)
            if (isPointOnSegment(d[0], a[0], b[0])) {
                return {{c, d}};
            }
            std::stringstream err;
            err.precision(std::numeric_limits<double>::digits10);
            err << ab << " x " << cd;
            OGS_FATAL(
                "The case of parallel line segments ({:s}) is not handled yet. "
                "Aborting.",
                err.str());
        }

        // check if d in (ab)
        if (isPointOnSegment(d[0], a[0], b[0])) {
            // check if a in (cd)
            if (isPointOnSegment(a[0], c[0], d[0])) {
                return {{a, d}};
            }
            // check if b==d
            if (MathLib::sqrDist2d(b, d) < eps) {
                return {{b}};
            }
            // check if b in (cd)
            if (isPointOnSegment(b[0], c[0], d[0])) {
                return {{b, d}};
            }
            // d in (ab), b not in (cd): check c in (ab)
            if (isPointOnSegment(c[0], a[0], b[0])) {
                return {{c, d}};
            }

            std::stringstream err;
            err.precision(std::numeric_limits<double>::digits10);
            err << ab << " x " << cd;
            OGS_FATAL(
                "The case of parallel line segments ({:s}) "
                "is not handled yet. Aborting.",
                err.str());
        }
        return std::vector<MathLib::Point3d>();
    }

    // precondition: points a, b, c are collinear
    // the function checks if the point c is onto the line segment (a,b)
    auto isCollinearPointOntoLineSegment = [](MathLib::Point3d const& a,
                                              MathLib::Point3d const& b,
                                              MathLib::Point3d const& c) {
        if (b[0] - a[0] != 0)
        {
            double const t = (c[0] - a[0]) / (b[0] - a[0]);
            return 0.0 <= t && t <= 1.0;
        }
        if (b[1] - a[1] != 0)
        {
            double const t = (c[1] - a[1]) / (b[1] - a[1]);
            return 0.0 <= t && t <= 1.0;
        }
        if (b[2] - a[2] != 0)
        {
            double const t = (c[2] - a[2]) / (b[2] - a[2]);
            return 0.0 <= t && t <= 1.0;
        }
        return false;
    };

    if (orient_abc == 0.0) {
        if (isCollinearPointOntoLineSegment(a, b, c))
        {
            return {{c}};
        }
        return std::vector<MathLib::Point3d>();
    }

    if (orient_abd == 0.0) {
        if (isCollinearPointOntoLineSegment(a, b, d))
        {
            return {{d}};
        }
        return std::vector<MathLib::Point3d>();
    }

    // check if the segment (ab) lies on the left or on the right of (cd)
    double const orient_cda(getOrientation(c, d, a));
    double const orient_cdb(getOrientation(c, d, b));
    if ((orient_cda > 0 && orient_cdb > 0) || (orient_cda < 0 && orient_cdb < 0)) {
        return std::vector<MathLib::Point3d>();
    }

    // at this point it is sure that there is an intersection and the system of
    // linear equations will be invertible
    // solve the two linear equations (b-a, c-d) (t, s)^T = (c-a) simultaneously
    Eigen::Matrix2d mat;
    mat(0,0) = b[0]-a[0];
    mat(0,1) = c[0]-d[0];
    mat(1,0) = b[1]-a[1];
    mat(1,1) = c[1]-d[1];
    Eigen::Vector2d rhs{c[0] - a[0], c[1] - a[1]};

    rhs = mat.partialPivLu().solve(rhs);
    if (0 <= rhs[1] && rhs[1] <= 1.0) {
        return { MathLib::Point3d{std::array<double,3>{{
                c[0]+rhs[1]*(d[0]-c[0]), c[1]+rhs[1]*(d[1]-c[1]),
                c[2]+rhs[1]*(d[2]-c[2])}} } };
    }
    return std::vector<MathLib::Point3d>();  // parameter s not in the valid
                                             // range
}

void sortSegments(
    MathLib::Point3d const& seg_beg_pnt,
    std::vector<GeoLib::LineSegment>& sub_segments)
{
    double const eps(std::numeric_limits<double>::epsilon());

    auto findNextSegment = [&eps](
                               MathLib::Point3d const& seg_beg_pnt,
                               std::vector<GeoLib::LineSegment>& sub_segments,
                               std::vector<GeoLib::LineSegment>::iterator&
                                   sub_seg_it) {
        if (sub_seg_it == sub_segments.end())
        {
            return;
        }
        // find appropriate segment for the given segment begin point
        auto act_beg_seg_it = std::find_if(
            sub_seg_it, sub_segments.end(),
            [&seg_beg_pnt, &eps](GeoLib::LineSegment const& seg)
            {
                return MathLib::sqrDist(seg_beg_pnt, seg.getBeginPoint()) < eps ||
                       MathLib::sqrDist(seg_beg_pnt, seg.getEndPoint()) < eps;
            });
        if (act_beg_seg_it == sub_segments.end())
        {
            return;
        }
        // if necessary correct orientation of segment, i.e. swap beg and end
        if (MathLib::sqrDist(seg_beg_pnt, act_beg_seg_it->getEndPoint()) <
            MathLib::sqrDist(seg_beg_pnt, act_beg_seg_it->getBeginPoint()))
        {
            std::swap(act_beg_seg_it->getBeginPoint(),
                      act_beg_seg_it->getEndPoint());
        }
        assert(sub_seg_it != sub_segments.end());
        // exchange segments within the container
        if (sub_seg_it != act_beg_seg_it)
        {
            std::swap(*sub_seg_it, *act_beg_seg_it);
        }
    };

    // find start segment
    auto seg_it = sub_segments.begin();
    findNextSegment(seg_beg_pnt, sub_segments, seg_it);

    while (seg_it != sub_segments.end())
    {
        MathLib::Point3d & new_seg_beg_pnt(seg_it->getEndPoint());
        seg_it++;
        if (seg_it != sub_segments.end())
        {
            findNextSegment(new_seg_beg_pnt, sub_segments, seg_it);
        }
    }
}

} // end namespace GeoLib
