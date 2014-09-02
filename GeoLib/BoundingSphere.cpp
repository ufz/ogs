/**
 * \file   Calculation of a minimum bounding sphere for a vector of points
 * \author Karsten Rink
 * \date   2014-07-11
 * \brief  Implementation of the BoundingSphere class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundingSphere.h"

#include <ctime>

#include "MathTools.h"
#include "AnalyticalGeometry.h"

namespace GeoLib {

BoundingSphere::BoundingSphere()
: _radius(-1), _center(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max())
{	
}

BoundingSphere::BoundingSphere(BoundingSphere const& sphere)
: _radius(sphere.getRadius()), _center(sphere.getCenter())
{
}

BoundingSphere::BoundingSphere(BoundingSphere const&& sphere)
: _radius(sphere.getRadius()), _center(sphere.getCenter())
{
}


BoundingSphere::BoundingSphere(GeoLib::Point const& p)
: _radius(std::numeric_limits<double>::epsilon()), _center(p)
{
}

BoundingSphere::BoundingSphere(GeoLib::Point const& p, double radius)
: _radius(radius), _center(p)
{
}

BoundingSphere::BoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q)
: _radius(std::numeric_limits<double>::epsilon()), _center(p)
{
    MathLib::Vector3 const a(p, q);

    if (a.getLength() > 0)
    {
        MathLib::Vector3 const o(0.5*a);
        _radius = o.getLength() + std::numeric_limits<double>::epsilon();
        _center = MathLib::Vector3(p) + o;
    }
}

BoundingSphere::BoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r)
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
        _radius = o.getLength() + std::numeric_limits<double>::epsilon();
        _center = MathLib::Vector3(p) + o;
    }
    else
    {
        BoundingSphere two_pnts_sphere;
        if (a.getLength() > b.getLength())
            two_pnts_sphere = BoundingSphere(p,r);
        else
            two_pnts_sphere = BoundingSphere(p,q);
        _radius = two_pnts_sphere.getRadius();
	    _center = two_pnts_sphere.getCenter();
    }
}

BoundingSphere::BoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r, GeoLib::Point const& s)
{
    MathLib::Vector3 const a(p, q);
    MathLib::Vector3 const b(p, r);
    MathLib::Vector3 const c(p, s);

    if (!GeoLib::isCoplanar(p, q, r, s))
    {
    	// det of matrix [a^T, b^T, c^T]^T
        double const denom = 2.0 * (a[0] * (b[1] * c[2] - c[1] * b[2])
                                  - b[0] * (a[1] * c[2] - c[1] * a[2])
                                  + c[0] * (a[1] * b[2] - b[1] * a[2]));
        MathLib::Vector3 const o = (scalarProduct(c,c) * crossProduct(a,b) 
                                  + scalarProduct(b,b) * crossProduct(c,a)
                                  + scalarProduct(a,a) * crossProduct(b,c)) 
                                  * (1.0 / denom);

        _radius = o.getLength() + std::numeric_limits<double>::epsilon();
        _center = MathLib::Vector3(p) + o;
    }
    else
    {
        BoundingSphere const pqr(p, q , r);
        BoundingSphere const pqs(p, q , s);
        BoundingSphere const prs(p, r , s);
        BoundingSphere const qrs(q, r , s);
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

BoundingSphere::BoundingSphere(std::vector<GeoLib::Point*> const& points)
: _center(0,0,0), _radius(-1)
{
	std::size_t const n_points (points.size());
	std::vector<GeoLib::Point*> sphere_points(points);
 
	BoundingSphere const bounding_sphere = recurseCalculation(sphere_points, 0, sphere_points.size(), 0);
	_center = bounding_sphere.getCenter();
	_radius = bounding_sphere.getRadius();
}

BoundingSphere BoundingSphere::recurseCalculation(std::vector<GeoLib::Point*> sphere_points, std::size_t start_idx, std::size_t length, std::size_t n_boundary_points)
{
    BoundingSphere sphere;
    switch(n_boundary_points)
    {
    case 0:
        sphere = BoundingSphere();
        break;
    case 1:
        sphere = BoundingSphere(*sphere_points[start_idx-1]);
        break;
    case 2:
        sphere = BoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2]);
        break;
    case 3:
        sphere = BoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2], *sphere_points[start_idx-3]);
        break;
    case 4:
        sphere = BoundingSphere(*sphere_points[start_idx-1], *sphere_points[start_idx-2], *sphere_points[start_idx-3], *sphere_points[start_idx-4]);
        return sphere;
    }

    for(std::size_t i=0; i<length; ++i)
    {
        // current point is located outside of sphere
        if (sphere.pointDistanceSquared(*sphere_points[start_idx+i]) > 0)
        {
            if (i>start_idx)
            {
                GeoLib::Point* tmp = sphere_points[start_idx+i];
                std::copy(sphere_points.begin() + start_idx, sphere_points.begin() + (start_idx + i), sphere_points.begin() + (start_idx + 1));
                sphere_points[start_idx] = tmp;
            }
            sphere = recurseCalculation(sphere_points, start_idx+1, i, n_boundary_points+1);
        }
    }
    return sphere;
}

double BoundingSphere::pointDistanceSquared(GeoLib::Point const& pnt) const
{
    return MathLib::sqrDist(_center.getCoords(), pnt.getCoords())-(_radius*_radius);
}

std::vector<GeoLib::Point*>* BoundingSphere::getRandomSpherePoints(std::size_t n_points) const
{
    std::vector<GeoLib::Point*> *pnts = new std::vector<GeoLib::Point*>;
    pnts->reserve(n_points);
    srand ( static_cast<unsigned>(time(NULL)) );

    for (std::size_t k(0); k<n_points; ++k) 
    {
        MathLib::Vector3 vec (0,0,0);
        double sum (0);
        for (unsigned i=0; i<3; ++i)
        {
            vec[i] = (double)rand()-(RAND_MAX/2.0);
            sum+=(vec[i]*vec[i]);
        }
        double const fac (_radius/sqrt(sum));
        pnts->push_back(new GeoLib::Point(_center[0]+vec[0]*fac, _center[1]+vec[1]*fac, _center[2]+vec[2]*fac));
    }
    return pnts;
}

}
