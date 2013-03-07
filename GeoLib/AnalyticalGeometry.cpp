/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Implementation of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AnalyticalGeometry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib> // for exit
#include <fstream>
#include <limits>
#include <list>

// BaseLib
#include "quicksort.h"

// GeoLib
#include "Polyline.h"
#include "Triangle.h"

// MathLib
#include "LinAlg/Solvers/GaussAlgorithm.h"
#include "MathTools.h"

namespace GeoLib
{
Orientation getOrientation(const double& p0_x, const double& p0_y, const double& p1_x,
                           const double& p1_y, const double& p2_x, const double& p2_y)
{
	double h1((p1_x - p0_x) * (p2_y - p0_y));
	double h2((p2_x - p0_x) * (p1_y - p0_y));

	double tol(sqrt( std::numeric_limits<double>::min()));
	if (fabs(h1 - h2) <= tol * std::max(fabs(h1), fabs(h2)))
		return COLLINEAR;
	if (h1 - h2 > 0.0)
		return CCW;

	return CW;
}

Orientation getOrientation(const GeoLib::Point* p0, const GeoLib::Point* p1,
                           const GeoLib::Point* p2)
{
	return getOrientation((*p0)[0], (*p0)[1], (*p1)[0], (*p1)[1], (*p2)[0], (*p2)[1]);
}

bool lineSegmentIntersect(const GeoLib::Point& a, const GeoLib::Point& b, const GeoLib::Point& c,
                          const GeoLib::Point& d, GeoLib::Point& s)
{
	MathLib::Matrix<double> mat(2, 2);
	mat(0, 0) = b[0] - a[0];
	mat(1, 0) = b[1] - a[1];
	mat(0,1) = c[0] - d[0];
	mat(1,1) = c[1] - d[1];

	// check if vectors are parallel
	double eps (sqrt(std::numeric_limits<double>::min()));
	if (fabs(mat(1,1)) < eps) {
		// vector (D-C) is parallel to x-axis
		if (fabs(mat(0,1)) < eps) {
			// vector (B-A) is parallel to x-axis
			return false;
		}
	} else {
		// vector (D-C) is not parallel to x-axis
		if (fabs(mat(0,1)) >= eps) {
			// vector (B-A) is not parallel to x-axis
			// \f$(B-A)\f$ and \f$(D-C)\f$ are parallel iff there exists
			// a constant \f$c\f$ such that \f$(B-A) = c (D-C)\f$
			if (fabs (mat(0,0) / mat(0,1) - mat(1,0) / mat(1,1)) < eps * fabs (mat(0,0) / mat(0,1)))
				return false;
		}
	}

	double *rhs (new double[2]);
	rhs[0] = c[0] - a[0];
	rhs[1] = c[1] - a[1];

	MathLib::GaussAlgorithm lu_solver (mat);
	lu_solver.execute (rhs);
	if (0 <= rhs[0] && rhs[0] <= 1.0 && 0 <= rhs[1] && rhs[1] <= 1.0) {
		s[0] = a[0] + rhs[0] * (b[0] - a[0]);
		s[1] = a[1] + rhs[0] * (b[1] - a[1]);
		s[2] = a[2] + rhs[0] * (b[2] - a[2]);
		// check z component
		double z0 (a[2] - d[2]), z1(rhs[0] * (b[2] - a[2]) + rhs[1] * (d[2] - c[2]));
		delete [] rhs;
		if (std::fabs (z0 - z1) < eps)
			return true;
		else
			return false;
	}
	else
		delete [] rhs;
	return false;
}

bool lineSegmentsIntersect(const GeoLib::Polyline* ply, size_t &idx0, size_t &idx1,
                           GeoLib::Point& intersection_pnt)
{
	size_t n_segs(ply->getNumberOfPoints() - 1);
	/**
	 * computing the intersections of all possible pairs of line segments of the given polyline
	 * as follows:
	 * let the segment \f$s_1 = (A,B)\f$ defined by \f$k\f$-th and \f$k+1\f$-st point
	 * of the polyline and segment \f$s_2 = (C,B)\f$ defined by \f$j\f$-th and
	 * \f$j+1\f$-st point of the polyline, \f$j>k+1\f$
	 */
	for (size_t k(0); k < n_segs - 2; k++) {
		for (size_t j(k + 2); j < n_segs; j++) {
			if (k != 0 || j < n_segs - 1) {
				if (lineSegmentIntersect(*(ply->getPoint(k)), *(ply->getPoint(k + 1)),
				                         *(ply->getPoint(j)), *(ply->getPoint(j + 1)),
				                         intersection_pnt)) {
					idx0 = k;
					idx1 = j;
					return true;
				}
			}
		}
	}
	return false;
}

bool isPointInTriangle(const double p[3], const double a[3], const double b[3], const double c[3])
{
	// criterion: p-b = u0 * (b - a) + u1 * (b - c); 0 <= u0, u1 <= 1, u0+u1 <= 1
	MathLib::Matrix<double> mat(2, 2);
	mat(0, 0) = a[0] - b[0];
	mat(0, 1) = c[0] - b[0];
	mat(1, 0) = a[1] - b[1];
	mat(1, 1) = c[1] - b[1];
	double rhs[2] = { p[0] - b[0], p[1] - b[1] };

	MathLib::GaussAlgorithm gauss(mat);
	gauss.execute(rhs);

	if (0 <= rhs[0] && rhs[0] <= 1 && 0 <= rhs[1] && rhs[1] <= 1 && rhs[0] + rhs[1] <= 1)
		return true;
	return false;
}

bool isPointInTriangle(const GeoLib::Point* p, const GeoLib::Point* a, const GeoLib::Point* b,
                       const GeoLib::Point* c)
{
	return isPointInTriangle(p->getCoords(), a->getCoords(), b->getCoords(), c->getCoords());
}

double getOrientedTriArea(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c)
{
	const double u[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
	const double v[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
	double w[3];
	MathLib::crossProd(u, v, w);
	return 0.5 * sqrt(MathLib::scpr<double, 3>(w, w));
}

bool isPointInTriangle(GeoLib::Point const& p, GeoLib::Point const& a, GeoLib::Point const& b,
                       GeoLib::Point const& c, double eps)
{
	const unsigned dim(3);
	MathLib::Matrix<double> m(dim, dim);
	for (unsigned i(0); i < dim; i++)
		m(i, 0) = b[i] - a[i];
	for (unsigned i(0); i < dim; i++)
		m(i, 1) = c[i] - a[i];
	for (unsigned i(0); i < dim; i++)
		m(i, 2) = p[i] - a[i];

	// point p is in the same plane as the triangle if and only if
	// the following determinate of the 3x3 matrix equals zero (up to an eps)
	double det3x3(m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2))
	              - m(1, 0) * (m(2, 1) * m(0, 2) - m(0, 1) * m(2, 2))
	              + m(2, 0) * (m(0, 1) * m(1, 2) - m(1, 1) * m(0, 2)));
	if (fabs(det3x3) > eps)
		return false;

	double total_area(getOrientedTriArea(a, b, c));
	double abp_area(getOrientedTriArea(a, b, p));
	double bcp_area(getOrientedTriArea(b, c, p));
	double cap_area(getOrientedTriArea(c, a, p));

	if (fabs(abp_area + bcp_area + cap_area - total_area) < eps)
		return true;
	return false;
}

// NewellPlane from book Real-Time Collision detection p. 494
void getNewellPlane(const std::vector<GeoLib::Point*>& pnts, MathLib::Vector &plane_normal, double& d)
{
	d = 0;
	MathLib::Vector centroid;
	size_t n_pnts(pnts.size());
	for (size_t i(n_pnts - 1), j(0); j < n_pnts; i = j, j++) {
		plane_normal[0] += ((*(pnts[i]))[1] - (*(pnts[j]))[1])
		                   * ((*(pnts[i]))[2] + (*(pnts[j]))[2]); // projection on yz
		plane_normal[1] += ((*(pnts[i]))[2] - (*(pnts[j]))[2])
		                   * ((*(pnts[i]))[0] + (*(pnts[j]))[0]); // projection on xz
		plane_normal[2] += ((*(pnts[i]))[0] - (*(pnts[j]))[0])
		                   * ((*(pnts[i]))[1] + (*(pnts[j]))[1]); // projection on xy

		centroid += *(pnts[j]);
	}

	plane_normal *= 1.0 / plane_normal.Length();
	d = centroid.Dot(plane_normal) / n_pnts;
}

void rotatePointsToXY(MathLib::Vector &plane_normal, std::vector<GeoLib::Point*> &pnts)
{
	double small_value(sqrt( std::numeric_limits<double>::min()));
	if (fabs(plane_normal[0]) < small_value && fabs(plane_normal[1]) < small_value)
		return;

	MathLib::Matrix<double> rot_mat(3, 3);
	computeRotationMatrixToXY(plane_normal, rot_mat);
	rotatePoints(rot_mat, pnts);

	double* tmp(rot_mat * plane_normal.getCoords());
	for (std::size_t j(0); j < 3; j++)
		plane_normal[j] = tmp[j];

	delete[] tmp;
}

void rotatePointsToXZ(MathLib::Vector &n, std::vector<GeoLib::Point*> &pnts)
{
	double small_value(sqrt( std::numeric_limits<double>::min()));
	if (fabs(n[0]) < small_value && fabs(n[1]) < small_value)
		return;

	// *** some frequently used terms ***
	// n_1^2 + n_2^2
	const double h0(n[0] * n[0] + n[1] * n[1]);
	// 1 / sqrt (n_1^2 + n_2^2)
	const double h1(1.0 / sqrt(h0));
	// 1 / sqrt (n_1^2 + n_2^2 + n_3^2)
	const double h2(1.0 / sqrt(h0 + n[2] * n[2]));

	MathLib::Matrix<double> rot_mat(3, 3);
	// calc rotation matrix
	rot_mat(0, 0) = n[1] * h1;
	rot_mat(0, 1) = -n[0] * h1;
	rot_mat(0, 2) = 0.0;
	rot_mat(1, 0) = n[0] * h2;
	rot_mat(1, 1) = n[1] * h2;
	rot_mat(1, 2) = n[2] * h2;
	rot_mat(2, 0) = n[0] * n[2] * h1 * h2;
	rot_mat(2, 1) = n[1] * n[2] * h1 * h2;
	rot_mat(2, 2) = -sqrt(h0) * h2;

	rotatePoints(rot_mat, pnts);

	double *tmp(rot_mat * n.getCoords());
	for (std::size_t j(0); j < 3; j++)
		n[j] = tmp[j];

	delete[] tmp;
}

void computeRotationMatrixToXY(MathLib::Vector const& plane_normal, MathLib::Matrix<double> & rot_mat)
{
	// *** some frequently used terms ***
	// sqrt (v_1^2 + v_2^2)
	double h0(sqrt(plane_normal[0] * plane_normal[0] + plane_normal[1]
	               * plane_normal[1]));
	// 1 / sqrt (v_1^2 + v_2^2)
	double h1(1 / h0);
	// 1 / sqrt (h0 + v_3^2)
	double h2(1.0 / sqrt(h0 + plane_normal[2] * plane_normal[2]));

	// calculate entries of rotation matrix
	rot_mat(0, 0) = plane_normal[2] * plane_normal[0] * h2 * h1;
	rot_mat(0, 1) = plane_normal[2] * plane_normal[1] * h2 * h1;
	rot_mat(0, 2) = -h0 * h2;
	rot_mat(1, 0) = -plane_normal[1] * h1;
	rot_mat(1, 1) = plane_normal[0] * h1;
	rot_mat(1, 2) = 0.0;
	rot_mat(2, 0) = plane_normal[0] * h2;
	rot_mat(2, 1) = plane_normal[1] * h2;
	rot_mat(2, 2) = plane_normal[2] * h2;
}

void rotatePoints(MathLib::Matrix<double> const& rot_mat, std::vector<GeoLib::Point*> &pnts)
{
	double* tmp (NULL);
	const std::size_t n_pnts(pnts.size());
	for (std::size_t k(0); k < n_pnts; k++) {
		tmp = rot_mat * pnts[k]->getCoords();
		for (std::size_t j(0); j < 3; j++)
			(*(pnts[k]))[j] = tmp[j];
		delete [] tmp;
	}
}
} // end namespace GeoLib
