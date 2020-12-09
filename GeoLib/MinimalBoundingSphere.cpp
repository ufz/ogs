/**
 * \file
 * \author Karsten Rink
 * \date   2014-07-11
 * \brief  Calculation of a minimum bounding sphere for a vector of points.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MinimalBoundingSphere.h"

#include <ctime>

#include "MathLib/Point3d.h"
#include "MathLib/GeometricBasics.h"
#include "MathLib/MathTools.h"

namespace GeoLib {
MinimalBoundingSphere::MinimalBoundingSphere() = default;

MinimalBoundingSphere::MinimalBoundingSphere(
    MathLib::Point3d const& p, double radius)
: _radius(radius), _center(p)
{
}

MinimalBoundingSphere::MinimalBoundingSphere(
    MathLib::Point3d const& p, MathLib::Point3d const& q)
: _radius(std::numeric_limits<double>::epsilon()), _center(p)
{
    auto const vp = Eigen::Map<Eigen::Vector3d const>(p.getCoords());
    auto const vq = Eigen::Map<Eigen::Vector3d const>(q.getCoords());
    Eigen::Vector3d const a = vq - vp;

    Eigen::Vector3d o = a / 2;
    _radius = o.norm() + std::numeric_limits<double>::epsilon();
    o += vp;
    _center = MathLib::Point3d{{o[0], o[1], o[2]}};
}

MinimalBoundingSphere::MinimalBoundingSphere(MathLib::Point3d const& p,
    MathLib::Point3d const& q, MathLib::Point3d const& r)
{
    auto const vp = Eigen::Map<Eigen::Vector3d const>(p.getCoords());
    auto const vq = Eigen::Map<Eigen::Vector3d const>(q.getCoords());
    auto const vr = Eigen::Map<Eigen::Vector3d const>(r.getCoords());
    Eigen::Vector3d const a = vr - vp;
    Eigen::Vector3d const b = vq - vp;
    Eigen::Vector3d const axb = a.cross(b);

    if (axb.squaredNorm() > 0)
    {
        double const denom = 2.0 * axb.dot(axb);
        Eigen::Vector3d o =
            (b.dot(b) * axb.cross(a) + a.dot(a) * b.cross(axb)) / denom;
        _radius = o.norm() + std::numeric_limits<double>::epsilon();
        o += vp;
        _center = MathLib::Point3d{{o[0], o[1], o[2]}};
    }
    else
    {
        MinimalBoundingSphere two_pnts_sphere;
        if (a.squaredNorm() > b.squaredNorm())
        {
            two_pnts_sphere = MinimalBoundingSphere(p,r);
        }
        else
        {
            two_pnts_sphere = MinimalBoundingSphere(p, q);
        }
        _radius = two_pnts_sphere.getRadius();
        _center = two_pnts_sphere.getCenter();
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(MathLib::Point3d const& p,
    MathLib::Point3d const& q,
    MathLib::Point3d const& r,
    MathLib::Point3d const& s)
{
    auto const vp = Eigen::Map<Eigen::Vector3d const>(p.getCoords());
    auto const vq = Eigen::Map<Eigen::Vector3d const>(q.getCoords());
    auto const vr = Eigen::Map<Eigen::Vector3d const>(r.getCoords());
    auto const vs = Eigen::Map<Eigen::Vector3d const>(s.getCoords());

    Eigen::Vector3d const va = vq - vp;
    Eigen::Vector3d const vb = vr - vp;
    Eigen::Vector3d const vc = vs - vp;

    if (!MathLib::isCoplanar(p, q, r, s))
    {
        double const denom = 2.0 * MathLib::scalarTriple(va, vb, vc);
        Eigen::Vector3d o =
            (vc.dot(vc) * va.cross(vb) + vb.dot(vb) * vc.cross(va) +
             va.dot(va) * vb.cross(vc)) /
            denom;

        _radius = o.norm() + std::numeric_limits<double>::epsilon();
        o += vp;
        _center = MathLib::Point3d({o[0], o[1], o[2]});
    }
    else
    {
        MinimalBoundingSphere const pqr(p, q , r);
        MinimalBoundingSphere const pqs(p, q , s);
        MinimalBoundingSphere const prs(p, r , s);
        MinimalBoundingSphere const qrs(q, r , s);
        _radius = pqr.getRadius();
        _center = pqr.getCenter();
        if (_radius < pqs.getRadius())
        {
            _radius = pqs.getRadius();
            _center = pqs.getCenter();
        }
        if (_radius < prs.getRadius())
        {
            _radius = prs.getRadius();
            _center = prs.getCenter();
        }
        if (_radius < qrs.getRadius())
        {
            _radius = qrs.getRadius();
            _center = qrs.getCenter();
        }
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(
    std::vector<MathLib::Point3d*> const& points)
    : _radius(-1), _center({0, 0, 0})
{
    const std::vector<MathLib::Point3d*>& sphere_points(points);
    MinimalBoundingSphere const bounding_sphere = recurseCalculation(sphere_points, 0, sphere_points.size(), 0);
    _center = bounding_sphere.getCenter();
    _radius = bounding_sphere.getRadius();
}

MinimalBoundingSphere
MinimalBoundingSphere::recurseCalculation(
    std::vector<MathLib::Point3d*> sphere_points,
    std::size_t start_idx,
    std::size_t length,
    std::size_t n_boundary_points)
{
    MinimalBoundingSphere sphere;
    switch(n_boundary_points)
    {
    case 0:
        sphere = MinimalBoundingSphere();
        break;
    case 1:
        sphere = MinimalBoundingSphere(*sphere_points[start_idx-1]);
        break;
    case 2:
        sphere = MinimalBoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2]);
        break;
    case 3:
        sphere = MinimalBoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2], *sphere_points[start_idx-3]);
        break;
    case 4:
        sphere = MinimalBoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2], *sphere_points[start_idx-3], *sphere_points[start_idx-4]);
        return sphere;
    }

    for(std::size_t i=0; i<length; ++i)
    {
        // current point is located outside of sphere
        if (sphere.pointDistanceSquared(*sphere_points[start_idx+i]) > 0)
        {
            if (i>start_idx)
            {
                using DiffType = std::vector<MathLib::Point3d*>::iterator::difference_type;
                std::vector<MathLib::Point3d*> const tmp_ps(
                    sphere_points.cbegin() + static_cast<DiffType>(start_idx),
                    sphere_points.cbegin() + static_cast<DiffType>(start_idx + i + 1));
                std::copy(tmp_ps.cbegin(), --tmp_ps.cend(),
                    sphere_points.begin() + static_cast<DiffType>(start_idx + 1));
                sphere_points[start_idx] = tmp_ps.back();
            }
            sphere = recurseCalculation(sphere_points, start_idx+1, i, n_boundary_points+1);
        }
    }
    return sphere;
}

double MinimalBoundingSphere::pointDistanceSquared(MathLib::Point3d const& pnt) const
{
    return MathLib::sqrDist(_center, pnt)-(_radius*_radius);
}

}  // namespace GeoLib
