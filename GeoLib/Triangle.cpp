/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-06
 * \brief  Implementation of the Triangle class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Triangle.h"

// GeoLib
#include "AnalyticalGeometry.h"

// MathLib
#include "LinAlg/Solvers/GaussAlgorithm.h"
#include "MathTools.h"
#include "LinAlg/Dense/DenseMatrix.h"
#include "Vector3.h"

namespace GeoLib {

Triangle::Triangle (std::vector<Point *> const &pnt_vec) :
	_pnts(pnt_vec), _initialized (false), _longest_edge (0.0)
{
	_pnt_ids[0] = std::numeric_limits<size_t>::max();
	_pnt_ids[1] = std::numeric_limits<size_t>::max();
	_pnt_ids[2] = std::numeric_limits<size_t>::max();
}

Triangle::Triangle (std::vector<Point *> const &pnt_vec, size_t pnt_a, size_t pnt_b, size_t pnt_c) :
	_pnts(pnt_vec), _initialized (true), _longest_edge (0.0)
{
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;
	_longest_edge = MathLib::sqrDist (*_pnts[_pnt_ids[0]], *_pnts[_pnt_ids[1]]);
	double tmp (MathLib::sqrDist (*_pnts[_pnt_ids[1]], *_pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = MathLib::sqrDist (*_pnts[_pnt_ids[0]], *_pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt (_longest_edge);
}

void Triangle::setTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert (pnt_a < _pnts.size() && pnt_b < _pnts.size() && pnt_c < _pnts.size());
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;

	_longest_edge = MathLib::sqrDist (*_pnts[_pnt_ids[0]], *_pnts[_pnt_ids[1]]);
	double tmp (MathLib::sqrDist (*_pnts[_pnt_ids[1]], *_pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = MathLib::sqrDist (*_pnts[_pnt_ids[0]], *_pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt (_longest_edge);
}

bool Triangle::containsPoint(Point const& q, double eps) const
{
	GeoLib::Point const& a(*(_pnts[_pnt_ids[0]]));
	GeoLib::Point const& b(*(_pnts[_pnt_ids[1]]));
	GeoLib::Point const& c(*(_pnts[_pnt_ids[2]]));
	return GeoLib::isPointInTriangle(q, a, b, c, eps);
}

bool Triangle::containsPoint2D (Point const& pnt) const
{
	GeoLib::Point const& a (*(_pnts[_pnt_ids[0]]));
	GeoLib::Point const& b (*(_pnts[_pnt_ids[1]]));
	GeoLib::Point const& c (*(_pnts[_pnt_ids[2]]));

	// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
	MathLib::DenseMatrix<double> mat (2,2);
	mat(0,0) = b[0] - a[0];
	mat(0,1) = c[0] - a[0];
	mat(1,0) = b[1] - a[1];
	mat(1,1) = c[1] - a[1];
	double y[2] = {pnt[0]-a[0], pnt[1]-a[1]};

	MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> gauss (mat);
	gauss.solve (y);

	const double delta (std::numeric_limits<double>::epsilon());
	const double upper (1+delta);

	// check if u0 and u1 fulfills the condition (with some delta)
	if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper && y[0] + y[1] <= upper) {
		return true;
	}
	return false;
}

void getPlaneCoefficients(Triangle const& tri, double c[3])
{
	GeoLib::Point const& p0 (*(tri.getPoint(0)));
	GeoLib::Point const& p1 (*(tri.getPoint(1)));
	GeoLib::Point const& p2 (*(tri.getPoint(2)));
	MathLib::DenseMatrix<double> mat (3,3);
	mat(0,0) = p0[0];
	mat(0,1) = p0[1];
	mat(0,2) = 1.0;
	mat(1,0) = p1[0];
	mat(1,1) = p1[1];
	mat(1,2) = 1.0;
	mat(2,0) = p2[0];
	mat(2,1) = p2[1];
	mat(2,2) = 1.0;
	c[0] = p0[2];
	c[1] = p1[2];
	c[2] = p2[2];

	MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> gauss (mat);
	gauss.solve (c);
}

} // end namespace GeoLib
