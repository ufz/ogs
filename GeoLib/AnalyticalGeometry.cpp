/**
 * \file
 * \author Thomas Fischer
 * \date   2010-03-17
 * \brief  Implementation of analytical geometry functions.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

	double tol(std::numeric_limits<double>::epsilon());
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

bool parallel(MathLib::Vector3 v, MathLib::Vector3 w)
{
	const double eps(std::numeric_limits<double>::epsilon());

	// check degenerated cases
	if (v.getLength() < eps)
		return false;

	if (w.getLength() < eps)
		return false;

	v.normalize();
	w.normalize();

	bool parallel(true);
	if (std::abs(v[0]-w[0]) > eps)
		parallel = false;
	if (std::abs(v[1]-w[1]) > eps)
		parallel = false;
	if (std::abs(v[2]-w[2]) > eps)
		parallel = false;

	if (! parallel) {
		parallel = true;
		// change sense of direction of v_normalised
		v *= -1.0;
		// check again
		if (std::abs(v[0]-w[0]) > eps)
			parallel = false;
		if (std::abs(v[1]-w[1]) > eps)
			parallel = false;
		if (std::abs(v[2]-w[2]) > eps)
			parallel = false;
	}

	return parallel;
}

bool lineSegmentIntersect(
	GeoLib::Point const& a,
	GeoLib::Point const& b,
	GeoLib::Point const& c,
	GeoLib::Point const& d,
	GeoLib::Point& s)
{
	if (!isCoplanar(a, b, c, d))
		return false;

	// handle special cases here to avoid computing intersection numerical
	if (MathLib::sqrDist(a, c) < std::numeric_limits<double>::epsilon() ||
		MathLib::sqrDist(a, d) < std::numeric_limits<double>::epsilon()) {
		s = a;
		return true;
	}
	if (MathLib::sqrDist(b, c) < std::numeric_limits<double>::epsilon() ||
		MathLib::sqrDist(b, d) < std::numeric_limits<double>::epsilon()) {
		s = b;
		return true;
	}

	// general case
	MathLib::Vector3 const v(a, b);
	MathLib::Vector3 const w(c, d);
	MathLib::Vector3 const qp(a, c);
	MathLib::Vector3 const pq(c, a);

	const double sqr_len_v(v.getSqrLength());
	const double sqr_len_w(w.getSqrLength());

	if (parallel(v,w)) {
		if (parallel(pq,v)) {
			const double sqr_dist_pq(pq.getSqrLength());
			if (sqr_dist_pq < sqr_len_v || sqr_dist_pq < sqr_len_w)
				return true;
		}
	}

	MathLib::DenseMatrix<double> mat(2,2);
	mat(0,0) = sqr_len_v;
	mat(0,1) = -1.0 * MathLib::scalarProduct(v,w);
	mat(1,1) = sqr_len_w;
	mat(1,0) = mat(0,1);

	double rhs[2] = {MathLib::scalarProduct(v,qp), MathLib::scalarProduct(w,pq)};

	MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> lu(mat);
	lu.solve (rhs);

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

	double const min_dist(sqrt(MathLib::sqrDist(p0, p1)));
	double const min_seg_len(std::min(sqrt(sqr_len_v), sqrt(sqr_len_w)));
	if (min_dist < min_seg_len * 1e-6) {
		s[0] = 0.5 * (p0[0] + p1[0]);
		s[1] = 0.5 * (p0[1] + p1[1]);
		s[2] = 0.5 * (p0[2] + p1[2]);
		return true;
	}

	return false;
}

bool lineSegmentsIntersect(const GeoLib::Polyline* ply,
                            size_t &idx0,
                            size_t &idx1,
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

bool isPointInTriangle(const GeoLib::Point* p, const GeoLib::Point* a, const GeoLib::Point* b,
                       const GeoLib::Point* c)
{
	return gaussPointInTriangle(*p, *a, *b, *c);
}

bool isPointInTriangle(GeoLib::Point const& p,
				GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
				double eps_pnt_out_of_plane, double eps_pnt_out_of_tri,
				GeoLib::TriangleTest algorithm)
{
	switch (algorithm)
	{
	case GeoLib::BARYCENTRIC:
		return barycentricPointInTriangle(p, a, b, c, eps_pnt_out_of_plane, eps_pnt_out_of_tri);
	default:
		return gaussPointInTriangle(p, a, b, c, eps_pnt_out_of_plane, eps_pnt_out_of_tri);
	}
	return false;
}

bool gaussPointInTriangle(GeoLib::Point const& q,
				GeoLib::Point const& a,
				GeoLib::Point const& b,
				GeoLib::Point const& c,
				double eps_pnt_out_of_plane, 
				double eps_pnt_out_of_tri)
{
	MathLib::Vector3 const v(a, b);
	MathLib::Vector3 const w(a, c);

	MathLib::DenseMatrix<double> mat (2,2);
	mat(0,0) = v.getSqrLength();
	mat(0,1) = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
	mat(1,0) = mat(0,1);
	mat(1,1) = w.getSqrLength();
	double y[2] = {
		v[0] * (q[0] - a[0]) + v[1] * (q[1] - a[1]) + v[2] * (q[2] - a[2]),
		w[0] * (q[0] - a[0]) + w[1] * (q[1] - a[1]) + w[2] * (q[2] - a[2])
	};

	MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> gauss(mat);
	gauss.solve(y);

	const double lower (eps_pnt_out_of_tri);
	const double upper (1 + lower);

	if (-lower <= y[0] && y[0] <= upper && -lower <= y[1] && y[1] <= upper && y[0] + y[1] <=
	    upper) {
		GeoLib::Point const q_projected(
			a[0] + y[0] * v[0] + y[1] * w[0],
			a[1] + y[0] * v[1] + y[1] * w[1],
			a[2] + y[0] * v[2] + y[1] * w[2]
		);
		if (MathLib::sqrDist(q, q_projected) < eps_pnt_out_of_plane)
			return true;
	}

	return false;
}

bool barycentricPointInTriangle(GeoLib::Point const& p,
				GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c,
				double eps_pnt_out_of_plane, double eps_pnt_out_of_tri)
{
	if (std::abs(orientation3d(p, a, b, c)) > eps_pnt_out_of_plane)
		return false;

	MathLib::Vector3 const pa (p,a); 
	MathLib::Vector3 const pb (p,b);
	MathLib::Vector3 const pc (p,c);
	double const area_x_2 (calcTriangleArea(a,b,c) * 2);
    
	double const alpha ((MathLib::crossProduct(pb,pc).getLength()) / area_x_2);
	if (alpha < -eps_pnt_out_of_tri || alpha > 1+eps_pnt_out_of_tri)
		return false;
	double const beta  ((MathLib::crossProduct(pc,pa).getLength()) / area_x_2);
	if (beta  < -eps_pnt_out_of_tri || beta  > 1+eps_pnt_out_of_tri)
		return false;
	double const gamma (1 - alpha - beta);
	if (gamma < -eps_pnt_out_of_tri || gamma > 1+eps_pnt_out_of_tri)
		return false;
	return true;
}

bool isPointInTetrahedron(GeoLib::Point const& p, GeoLib::Point const& a, GeoLib::Point const& b, 
                          GeoLib::Point const& c, GeoLib::Point const& d, double eps)
{
    double const d0 (orientation3d(d,a,b,c));
    // if tetrahedron is not coplanar
    if (std::abs(d0) > std::numeric_limits<double>::epsilon())
    {
        bool const d0_sign (d0>0);
        // if p is on the same side of bcd as a
        double const d1 (orientation3d(d, p, b, c));
        if (!(d0_sign == (d1>=0) || std::abs(d1) < eps))
            return false;
        // if p is on the same side of acd as b
        double const d2 (orientation3d(d, a, p, c));
        if (!(d0_sign == (d2>=0) || std::abs(d2) < eps))
            return false;
        // if p is on the same side of abd as c
        double const d3 (orientation3d(d, a, b, p));
        if (!(d0_sign == (d3>=0) || std::abs(d3) < eps))
            return false;
        // if p is on the same side of abc as d
        double const d4 (orientation3d(p, a, b, c));
        if (!(d0_sign == (d4>=0) || std::abs(d4) < eps))
            return false;
        return true;
    }
    return false;
}

double calcTriangleArea(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c)
{
	MathLib::Vector3 const u(a,c);
	MathLib::Vector3 const v(a,b);
	MathLib::Vector3 const w(MathLib::crossProduct(u, v));
	return 0.5 * w.getLength();
}

double calcTetrahedronVolume(const double* x1, const double* x2, const double* x3, const double* x4)
{
	const MathLib::Vector3 ab(x1, x2);
	const MathLib::Vector3 ac(x1, x3);
	const MathLib::Vector3 ad(x1, x4);
	return std::abs(GeoLib::scalarTriple(ac, ad, ab)) / 6.0;
}

// NewellPlane from book Real-Time Collision detection p. 494
void getNewellPlane(const std::vector<GeoLib::Point*>& pnts, MathLib::Vector3 &plane_normal, double& d)
{
	d = 0;
	MathLib::Vector3 centroid;
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

	plane_normal *= 1.0 / plane_normal.getLength();
	d = MathLib::scalarProduct(centroid, plane_normal) / n_pnts;
}

void computeRotationMatrixToXY(MathLib::Vector3 const& plane_normal, MathLib::DenseMatrix<double> & rot_mat)
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

void computeRotationMatrixToXZ(MathLib::Vector3 const& plane_normal, MathLib::DenseMatrix<double> & rot_mat)
{
	// *** some frequently used terms ***
	// n_1^2 + n_2^2
	const double h0(plane_normal[0] * plane_normal[0] + plane_normal[1] * plane_normal[1]);
	// 1 / sqrt (n_1^2 + n_2^2)
	const double h1(1.0 / sqrt(h0));
	// 1 / sqrt (n_1^2 + n_2^2 + n_3^2)
	const double h2(1.0 / sqrt(h0 + plane_normal[2] * plane_normal[2]));

	// calc rotation matrix
	rot_mat(0, 0) = plane_normal[1] * h1;
	rot_mat(0, 1) = -plane_normal[0] * h1;
	rot_mat(0, 2) = 0.0;
	rot_mat(1, 0) = plane_normal[0] * h2;
	rot_mat(1, 1) = plane_normal[1] * h2;
	rot_mat(1, 2) = plane_normal[2] * h2;
	rot_mat(2, 0) = plane_normal[0] * plane_normal[2] * h1 * h2;
	rot_mat(2, 1) = plane_normal[1] * plane_normal[2] * h1 * h2;
	rot_mat(2, 2) = -sqrt(h0) * h2;
}

void rotatePoints(MathLib::DenseMatrix<double> const& rot_mat, std::vector<GeoLib::Point*> &pnts)
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

void rotatePointsToXY(std::vector<GeoLib::Point*> &pnts)
{
	assert(pnts.size()>2);
	// calculate supporting plane
	MathLib::Vector3 plane_normal;
	double d;
	// compute the plane normal
	GeoLib::getNewellPlane(pnts, plane_normal, d);

	const double tol (std::numeric_limits<double>::epsilon());
	if (std::abs(plane_normal[0]) > tol || std::abs(plane_normal[1]) > tol) {
		// rotate copied points into x-y-plane
		MathLib::DenseMatrix<double> rot_mat(3, 3);
		computeRotationMatrixToXY(plane_normal, rot_mat);
		rotatePoints(rot_mat, pnts);
	}

	for (std::size_t k(0); k<pnts.size(); k++)
		(*(pnts[k]))[2] = 0.0; // should be -= d but there are numerical errors
}

void rotatePointsToXZ(std::vector<GeoLib::Point*> &pnts)
{
	assert(pnts.size()>2);
	// calculate supporting plane
	MathLib::Vector3 plane_normal;
	double d;
	// compute the plane normal
	GeoLib::getNewellPlane(pnts, plane_normal, d);

	const double tol (std::numeric_limits<double>::epsilon());
	if (std::abs(plane_normal[0]) > tol || std::abs(plane_normal[1]) > tol) {
		// rotate copied points into x-z-plane
		MathLib::DenseMatrix<double> rot_mat(3, 3);
		computeRotationMatrixToXZ(plane_normal, rot_mat);
		rotatePoints(rot_mat, pnts);
	}

	for (std::size_t k(0); k<pnts.size(); k++)
		(*(pnts[k]))[1] = 0.0; // should be -= d but there are numerical errors
}

GeoLib::Point* triangleLineIntersection(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c, GeoLib::Point const& p, GeoLib::Point const& q)
{
	const MathLib::Vector3 pq(p, q);
	const MathLib::Vector3 pa(p, a);
	const MathLib::Vector3 pb(p, b);
	const MathLib::Vector3 pc(p, c);
	
	double u (scalarTriple(pq, pc, pb));
	if (u<0) return nullptr;
	double v (scalarTriple(pq, pa, pc));
	if (v<0) return nullptr;
	double w (scalarTriple(pq, pb, pa));
	if (w<0) return nullptr;
	
	const double denom (1.0/(u+v+w));
	u*=denom;
	v*=denom;
	w*=denom;
	return new GeoLib::Point(u*a[0]+v*b[0]+w*c[0],u*a[1]+v*b[1]+w*c[1],u*a[2]+v*b[2]+w*c[2]);
}

double scalarTriple(MathLib::Vector3 const& u, MathLib::Vector3 const& v, MathLib::Vector3 const& w)
{
	MathLib::Vector3 const cross(MathLib::crossProduct(u, v));
	return MathLib::scalarProduct(cross,w);
}

double orientation3d(GeoLib::Point const& p,
                     GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c)
{
    MathLib::Vector3 const ap (a, p);
    MathLib::Vector3 const bp (b, p);
    MathLib::Vector3 const cp (c, p);
    return scalarTriple(bp,cp,ap);
}

bool dividedByPlane(const GeoLib::Point& a, const GeoLib::Point& b, const GeoLib::Point& c, const GeoLib::Point& d)
{
	for (unsigned x=0; x<3; ++x)
	{
		const unsigned y=(x+1)%3;
		const double abc = (b[x] - a[x])*(c[y] - a[y]) - (b[y] - a[y])*(c[x] - a[x]);
		const double abd = (b[x] - a[x])*(d[y] - a[y]) - (b[y] - a[y])*(d[x] - a[x]);

		if ((abc>0 && abd<0) || (abc<0 && abd>0))
			return true;		
	}
	return false;
}

bool isCoplanar(const GeoLib::Point& a, const GeoLib::Point& b, const GeoLib::Point& c, const GeoLib::Point& d)
{
	const MathLib::Vector3 ab(a,b);
	const MathLib::Vector3 ac(a,c);
	const MathLib::Vector3 ad(a,d);

	if (ab.getSqrLength() < pow(std::numeric_limits<double>::epsilon(),2) ||
		ac.getSqrLength() < pow(std::numeric_limits<double>::epsilon(),2) ||
		ad.getSqrLength() < pow(std::numeric_limits<double>::epsilon(),2)) {
		return true;
	}

	// In exact arithmetic <ac*ad^T, ab> should be zero
	// if all four points are coplanar.
	const double sqr_scalar_triple(pow(MathLib::scalarProduct(MathLib::crossProduct(ac,ad), ab),2));
	// Due to evaluating the above numerically some cancellation or rounding
	// can occur. For this reason a normalisation factor is introduced.
	const double normalisation_factor =
		(ab.getSqrLength() * ac.getSqrLength() * ad.getSqrLength());

	// tolerance 1e-11 is choosen such that
	// a = (0,0,0), b=(1,0,0), c=(0,1,0) and d=(1,1,1e-6) are considered as coplanar
	// a = (0,0,0), b=(1,0,0), c=(0,1,0) and d=(1,1,1e-5) are considered as not coplanar
	return (sqr_scalar_triple/normalisation_factor < 1e-11);
}

void computeAndInsertAllIntersectionPoints(GeoLib::PointVec &pnt_vec,
	std::vector<GeoLib::Polyline*> & plys)
{
	for (auto it0(plys.begin()); it0 != plys.end(); ++it0) {
		auto it1(it0);
		++it1;
		for (; it1 != plys.end(); ++it1) {
			GeoLib::Point s;
			for (std::size_t i(0); i<(*it0)->getNumberOfPoints()-1; i++) {
				for (std::size_t j(0); j<(*it1)->getNumberOfPoints()-1; j++) {
					if (lineSegmentIntersect(*(*it0)->getPoint(i), *(*it0)->getPoint(i+1),
						*(*it1)->getPoint(j), *(*it1)->getPoint(j+1), s)) {
						std::size_t const id(pnt_vec.push_back(new GeoLib::Point(s)));
						(*it0)->insertPoint(i+1, id);
						(*it1)->insertPoint(j+1, id);
					}
				}
			}
		}
	}
}


} // end namespace GeoLib
