/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Triangle.h"

#include "Point.h"
#include "AnalyticalGeometry.h"

#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/GeometricBasics.h"

namespace GeoLib {

Triangle::Triangle (std::vector<Point *> const &pnt_vec,
    std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c) :
    _pnts(pnt_vec), _pnt_ids( {{pnt_a, pnt_b, pnt_c}} )
{
    assert(!_pnts.empty());
    assert (pnt_a < _pnts.size() && pnt_b < _pnts.size() && pnt_c < _pnts.size());
}

void Triangle::setTriangle (std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c)
{
    assert (pnt_a < _pnts.size() && pnt_b < _pnts.size() && pnt_c < _pnts.size());
    _pnt_ids[0] = pnt_a;
    _pnt_ids[1] = pnt_b;
    _pnt_ids[2] = pnt_c;
}

bool Triangle::containsPoint(MathLib::Point3d const& q, double eps) const
{
    GeoLib::Point const& a(*(_pnts[_pnt_ids[0]]));
    GeoLib::Point const& b(*(_pnts[_pnt_ids[1]]));
    GeoLib::Point const& c(*(_pnts[_pnt_ids[2]]));
    return MathLib::isPointInTriangle(q, a, b, c, eps);
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

    MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> gauss;
    gauss.solve(mat, y, true);

    const double delta (std::numeric_limits<double>::epsilon());
    const double upper (1+delta);

    // check if u0 and u1 fulfills the condition (with some delta)
    if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper &&
        y[0] + y[1] <= upper)
    {
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

    MathLib::GaussAlgorithm<MathLib::DenseMatrix<double>, double*> gauss;
    gauss.solve (mat, c);
}

} // end namespace GeoLib
