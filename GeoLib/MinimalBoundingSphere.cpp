/**
 * \file   Calculation of a minimum bounding sphere for a vector of points
 * \author Karsten Rink
 * \date   2014-07-11
 * \brief  Implementation of the MinimalBoundingSphere class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MinimalBoundingSphere.h"

#include <ctime>

#include "MathTools.h"
#include "AnalyticalGeometry.h"

namespace GeoLib {

MinimalBoundingSphere::MinimalBoundingSphere()
: _radius(-1), _center(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max())
{	
}

MinimalBoundingSphere::MinimalBoundingSphere(GeoLib::Point const& p, double radius)
: _radius(radius), _center(p)
{
}

MinimalBoundingSphere::MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q)
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

MinimalBoundingSphere::MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r)
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
        MinimalBoundingSphere two_pnts_sphere;
        if (a.getLength() > b.getLength())
            two_pnts_sphere = MinimalBoundingSphere(p,r);
        else
            two_pnts_sphere = MinimalBoundingSphere(p,q);
        _radius = two_pnts_sphere.getRadius();
	    _center = two_pnts_sphere.getCenter();
    }
}

MinimalBoundingSphere::MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r, GeoLib::Point const& s)
{
    MathLib::Vector3 const a(p, q);
    MathLib::Vector3 const b(p, r);
    MathLib::Vector3 const c(p, s);

    if (!GeoLib::isCoplanar(p, q, r, s))
    {
        double const denom = 2.0 * GeoLib::scalarTriple(a,b,c);
        MathLib::Vector3 const o = (scalarProduct(c,c) * crossProduct(a,b) 
                                  + scalarProduct(b,b) * crossProduct(c,a)
                                  + scalarProduct(a,a) * crossProduct(b,c)) 
                                  * (1.0 / denom);

        _radius = o.getLength() + std::numeric_limits<double>::epsilon();
        _center = MathLib::Vector3(p) + o;
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

MinimalBoundingSphere::MinimalBoundingSphere(std::vector<GeoLib::Point*> const& points)
: _radius(-1), _center(0,0,0)
{
	std::vector<GeoLib::Point*> sphere_points(points);
 	MinimalBoundingSphere const bounding_sphere = recurseCalculation(sphere_points, 0, sphere_points.size(), 0);
	_center = bounding_sphere.getCenter();
	_radius = bounding_sphere.getRadius();
}

MinimalBoundingSphere MinimalBoundingSphere::recurseCalculation(std::vector<GeoLib::Point*> sphere_points, std::size_t start_idx, std::size_t length, std::size_t n_boundary_points)
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
                GeoLib::Point* tmp = sphere_points[start_idx+i];
                std::copy(sphere_points.begin() + start_idx, sphere_points.begin() + (start_idx + i), sphere_points.begin() + (start_idx + 1));
                sphere_points[start_idx] = tmp;
            }
            sphere = recurseCalculation(sphere_points, start_idx+1, i, n_boundary_points+1);
        }
    }
    return sphere;
}

double MinimalBoundingSphere::pointDistanceSquared(GeoLib::Point const& pnt) const
{
    return MathLib::sqrDist(_center.getCoords(), pnt.getCoords())-(_radius*_radius);
}

std::vector<GeoLib::Point*>* MinimalBoundingSphere::getRandomSpherePoints(std::size_t n_points) const
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
