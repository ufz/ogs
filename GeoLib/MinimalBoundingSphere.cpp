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
#include "MathLib/Vector3.h"

namespace GeoLib {
MinimalBoundingSphere::MinimalBoundingSphere() = default;

MinimalBoundingSphere::MinimalBoundingSphere(
    MathLib::Point3d const& p, double radius)
: radius_(radius), center_(p)
{
}

MinimalBoundingSphere::MinimalBoundingSphere(
    MathLib::Point3d const& p, MathLib::Point3d const& q)
: radius_(std::numeric_limits<double>::epsilon()), center_(p)
{
    MathLib::Vector3 const a(p, q);

    if (a.getLength() > 0)
    {
        MathLib::Vector3 const o(0.5*a);
        radius_ = o.getLength() + std::numeric_limits<double>::epsilon();
        center_ = MathLib::Vector3(p) + o;
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(MathLib::Point3d const& p,
    MathLib::Point3d const& q, MathLib::Point3d const& r)
{
    MathLib::Vector3 const a(p,r);
    MathLib::Vector3 const b(p,q);

    MathLib::Vector3 const cross_ab(crossProduct(a,b));

    if (cross_ab.getLength() > 0)
    {
        double const denom = 2.0 * scalarProduct(cross_ab,cross_ab);
        MathLib::Vector3 const o = (scalarProduct(b,b) * crossProduct(cross_ab, a)
                                   + scalarProduct(a,a) * crossProduct(b, cross_ab))
                                  * (1.0 / denom);
        radius_ = o.getLength() + std::numeric_limits<double>::epsilon();
        center_ = MathLib::Vector3(p) + o;
    }
    else
    {
        MinimalBoundingSphere two_pnts_sphere;
        if (a.getLength() > b.getLength())
        {
            two_pnts_sphere = MinimalBoundingSphere(p,r);
        }
        else
        {
            two_pnts_sphere = MinimalBoundingSphere(p, q);
        }
        radius_ = two_pnts_sphere.getRadius();
        center_ = MathLib::Vector3(two_pnts_sphere.getCenter());
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(MathLib::Point3d const& p,
    MathLib::Point3d const& q,
    MathLib::Point3d const& r,
    MathLib::Point3d const& s)
{
    MathLib::Vector3 const a(p, q);
    MathLib::Vector3 const b(p, r);
    MathLib::Vector3 const c(p, s);

    if (!MathLib::isCoplanar(p, q, r, s))
    {
        double const denom = 2.0 * MathLib::scalarTriple(a,b,c);
        MathLib::Vector3 const o = (scalarProduct(c,c) * crossProduct(a,b)
                                  + scalarProduct(b,b) * crossProduct(c,a)
                                  + scalarProduct(a,a) * crossProduct(b,c))
                                  * (1.0 / denom);

        radius_ = o.getLength() + std::numeric_limits<double>::epsilon();
        center_ = MathLib::Vector3(p) + o;
    }
    else
    {
        MinimalBoundingSphere const pqr(p, q , r);
        MinimalBoundingSphere const pqs(p, q , s);
        MinimalBoundingSphere const prs(p, r , s);
        MinimalBoundingSphere const qrs(q, r , s);
        radius_ = pqr.getRadius();
        center_ = MathLib::Vector3(pqr.getCenter());
        if (radius_ < pqs.getRadius())
        {
            radius_ = pqs.getRadius();
            center_ = MathLib::Vector3(pqs.getCenter());
        }
        if (radius_ < prs.getRadius())
        {
            radius_ = prs.getRadius();
            center_ = MathLib::Vector3(prs.getCenter());
        }
        if (radius_ < qrs.getRadius())
        {
            radius_ = qrs.getRadius();
            center_ = MathLib::Vector3(qrs.getCenter());
        }
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(
    std::vector<MathLib::Point3d*> const& points)
: radius_(-1), center_(0,0,0)
{
    const std::vector<MathLib::Point3d*>& sphere_points(points);
    MinimalBoundingSphere const bounding_sphere = recurseCalculation(sphere_points, 0, sphere_points.size(), 0);
    center_ = MathLib::Vector3(bounding_sphere.getCenter());
    radius_ = bounding_sphere.getRadius();
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
    return MathLib::sqrDist(center_.getCoords(), pnt.getCoords())-(radius_*radius_);
}

std::vector<MathLib::Point3d*>* MinimalBoundingSphere::getRandomSpherePoints(std::size_t n_points) const
{
    auto* pnts = new std::vector<MathLib::Point3d*>;
    pnts->reserve(n_points);
    srand ( static_cast<unsigned>(time(nullptr)) );

    for (std::size_t k(0); k<n_points; ++k)
    {
        MathLib::Vector3 vec (0,0,0);
        double sum (0);
        for (unsigned i=0; i<3; ++i)
        {
            vec[i] = static_cast<double>(rand())-(RAND_MAX/2.0);
            sum+=(vec[i]*vec[i]);
        }
        double const fac (radius_/sqrt(sum));
        pnts->push_back(new MathLib::Point3d(center_+fac * vec));
    }
    return pnts;
}

}  // namespace GeoLib
