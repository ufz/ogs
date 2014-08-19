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

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MathTools.h"

namespace GeoLib {

BoundingSphere::BoundingSphere()
: _center(0,0,0), _radius(-1)
{	
}

BoundingSphere::BoundingSphere(const BoundingSphere &sphere)
: _center(sphere.getCenter()), _radius(sphere.getRadius())
{
}

BoundingSphere::BoundingSphere(const GeoLib::Point &p)
: _center(p), _radius(std::numeric_limits<double>::epsilon())
{
}

BoundingSphere::BoundingSphere(const GeoLib::Point &p, double radius)
: _center(p), _radius(radius)
{
}

BoundingSphere::BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q)
{
	const MathLib::Vector3 a(p, q);
	const MathLib::Vector3 o(0.5*a);
	_radius = o.getLength() + std::numeric_limits<double>::epsilon();
	_center = MathLib::Vector3(p) + o;
}

BoundingSphere::BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q,  const GeoLib::Point &r)
{
	const MathLib::Vector3 a(p,r);
	const MathLib::Vector3 b(p,q);

	const MathLib::Vector3 cross_ab(crossProduct(a,b));
	const double denom = 2.0 * scalarProduct(cross_ab,cross_ab);
	const MathLib::Vector3 o = (scalarProduct(b,b) * crossProduct(cross_ab, a) 
	 	                      + scalarProduct(a,a) * crossProduct(b, cross_ab))
		                      * (1.0 / denom);
	_radius = o.getLength() + std::numeric_limits<double>::epsilon();
	_center = MathLib::Vector3(p) + o;
}

BoundingSphere::BoundingSphere(const GeoLib::Point &p, const GeoLib::Point &q, const GeoLib::Point &r, const GeoLib::Point &s)
{
	const MathLib::Vector3 a(p, q);
	const MathLib::Vector3 b(p, r);
	const MathLib::Vector3 c(p, s);

	// det of matrix [a^T, b^T, c^T]^T
	const double denom = 2.0 * (a[0] * (b[1] * c[2] - c[1] * b[2])
							  - b[0] * (a[1] * c[2] - c[1] * a[2])
							  + c[0] * (a[1] * b[2] - b[1] * a[2]));
	const MathLib::Vector3 o = (scalarProduct(c,c) * crossProduct(a,b) 
							  + scalarProduct(b,b) * crossProduct(c,a)
							  + scalarProduct(a,a) * crossProduct(b,c)) 
							* (1.0 / denom);

	_radius = o.getLength() + std::numeric_limits<double>::epsilon();
	_center = MathLib::Vector3(p) + o;
}

BoundingSphere::BoundingSphere(const std::vector<GeoLib::Point*> &points)
: _center(0,0,0), _radius(-1)
{
	const std::size_t n_points (points.size());
	GeoLib::Point **sphere_points = new GeoLib::Point*[n_points];
	for(unsigned int i = 0; i < n_points; i++)
		sphere_points[i] = points[i];

	const BoundingSphere bounding_sphere = recurseCalculation(sphere_points, n_points, 0);
	delete[] sphere_points;
	
	this->_center = bounding_sphere.getCenter();
	this->_radius = bounding_sphere.getRadius();
}

BoundingSphere BoundingSphere::recurseCalculation(GeoLib::Point* sphere_points[], std::size_t n_points, std::size_t boundary_points)
{
	BoundingSphere sphere;
	switch(boundary_points)
	{
	case 0:
		sphere = BoundingSphere();
		break;
	case 1:
		sphere = BoundingSphere(*sphere_points[-1]);
		break;
	case 2:
		sphere = BoundingSphere(*sphere_points[-1], *sphere_points[-2]);
		break;
	case 3:
		sphere = BoundingSphere(*sphere_points[-1], *sphere_points[-2], *sphere_points[-3]);
		break;
	case 4:
		sphere = BoundingSphere(*sphere_points[-1], *sphere_points[-2], *sphere_points[-3], *sphere_points[-4]);
		return sphere;
	}

	for(std::size_t i=0; i<n_points; ++i)
	{
		if(sphere.sqrPointDist(*sphere_points[i]) > 0)
		{
			for(unsigned int j = i; j > 0; j--)
			{
				GeoLib::Point* tmp = sphere_points[j];
				sphere_points[j] = sphere_points[j - 1];
				sphere_points[j - 1] = tmp;
			}
		}
		sphere = recurseCalculation(sphere_points+1, i, boundary_points+1);
	}
	return sphere;
}

double BoundingSphere::sqrPointDist(const GeoLib::Point pnt) const
{
	return MathLib::sqrDist(_center.getCoords(), pnt.getCoords())-(_radius*_radius);
}

std::vector<GeoLib::Point*>* BoundingSphere::getSpherePoints(std::size_t n_points) const
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
		double fac (this->_radius/sqrt(sum));
		pnts->push_back(new GeoLib::Point(_center[0]+vec[0]*fac, _center[1]+vec[1]*fac, _center[2]+vec[2]*fac));
	}
	return pnts;
}

}
