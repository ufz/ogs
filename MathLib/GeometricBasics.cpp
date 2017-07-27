/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <logog/include/logog.hpp>

#include <Eigen/Dense>

#include "Point3d.h"
#include "Vector3.h"

#include "GeometricBasics.h"

namespace MathLib
{
double orientation3d(MathLib::Point3d const& p,
                     MathLib::Point3d const& a,
                     MathLib::Point3d const& b,
                     MathLib::Point3d const& c)
{
    MathLib::Vector3 const ap (a, p);
    MathLib::Vector3 const bp (b, p);
    MathLib::Vector3 const cp (c, p);
    return MathLib::scalarTriple(bp,cp,ap);
}

double calcTetrahedronVolume(MathLib::Point3d const& a,
                             MathLib::Point3d const& b,
                             MathLib::Point3d const& c,
                             MathLib::Point3d const& d)
{
    const MathLib::Vector3 ab(a, b);
    const MathLib::Vector3 ac(a, c);
    const MathLib::Vector3 ad(a, d);
    return std::abs(MathLib::scalarTriple(ac, ad, ab)) / 6.0;
}

double calcTriangleArea(MathLib::Point3d const& a, MathLib::Point3d const& b,
                        MathLib::Point3d const& c)
{
    MathLib::Vector3 const u(a, c);
    MathLib::Vector3 const v(a, b);
    MathLib::Vector3 const w(MathLib::crossProduct(u, v));
    return 0.5 * w.getLength();
}

bool isPointInTetrahedron(MathLib::Point3d const& p, MathLib::Point3d const& a,
                          MathLib::Point3d const& b, MathLib::Point3d const& c,
                          MathLib::Point3d const& d, double eps)
{
    double const d0 (MathLib::orientation3d(d,a,b,c));
    // if tetrahedron is not coplanar
    if (std::abs(d0) > std::numeric_limits<double>::epsilon())
    {
        bool const d0_sign (d0>0);
        // if p is on the same side of bcd as a
        double const d1 (MathLib::orientation3d(d, p, b, c));
        if (!(d0_sign == (d1>=0) || std::abs(d1) < eps))
            return false;
        // if p is on the same side of acd as b
        double const d2 (MathLib::orientation3d(d, a, p, c));
        if (!(d0_sign == (d2>=0) || std::abs(d2) < eps))
            return false;
        // if p is on the same side of abd as c
        double const d3 (MathLib::orientation3d(d, a, b, p));
        if (!(d0_sign == (d3>=0) || std::abs(d3) < eps))
            return false;
        // if p is on the same side of abc as d
        double const d4 (MathLib::orientation3d(p, a, b, c));
        return d0_sign == (d4 >= 0) || std::abs(d4) < eps;
    }
    return false;
}

bool isPointInTriangle(MathLib::Point3d const& p,
                       MathLib::Point3d const& a,
                       MathLib::Point3d const& b,
                       MathLib::Point3d const& c,
                       double eps_pnt_out_of_plane,
                       double eps_pnt_out_of_tri,
                       MathLib::TriangleTest algorithm)
{
    switch (algorithm)
    {
        case MathLib::GAUSS:
            return gaussPointInTriangle(p, a, b, c, eps_pnt_out_of_plane,
                                        eps_pnt_out_of_tri);
        case MathLib::BARYCENTRIC:
            return barycentricPointInTriangle(p, a, b, c, eps_pnt_out_of_plane,
                                              eps_pnt_out_of_tri);
        default:
            ERR("Selected algorithm for point in triangle testing not found, "
                "falling back on default.");
    }
    return gaussPointInTriangle(p, a, b, c, eps_pnt_out_of_plane,
                                eps_pnt_out_of_tri);
}

bool gaussPointInTriangle(MathLib::Point3d const& q,
                          MathLib::Point3d const& a,
                          MathLib::Point3d const& b,
                          MathLib::Point3d const& c,
                          double eps_pnt_out_of_plane,
                          double eps_pnt_out_of_tri)
{
    MathLib::Vector3 const v(a, b);
    MathLib::Vector3 const w(a, c);

    Eigen::Matrix2d mat;
    mat(0, 0) = v.getSqrLength();
    mat(0, 1) = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    mat(1, 0) = mat(0, 1);
    mat(1, 1) = w.getSqrLength();
    Eigen::Vector2d y;
    y << v[0] * (q[0] - a[0]) + v[1] * (q[1] - a[1]) + v[2] * (q[2] - a[2]),
        w[0] * (q[0] - a[0]) + w[1] * (q[1] - a[1]) + w[2] * (q[2] - a[2]);

    y = mat.partialPivLu().solve(y);

    const double lower(eps_pnt_out_of_tri);
    const double upper(1 + lower);

    if (-lower <= y[0] && y[0] <= upper && -lower <= y[1] && y[1] <= upper &&
        y[0] + y[1] <= upper)
    {
        MathLib::Point3d const q_projected(std::array<double, 3>{
            {a[0] + y[0] * v[0] + y[1] * w[0], a[1] + y[0] * v[1] + y[1] * w[1],
             a[2] + y[0] * v[2] + y[1] * w[2]}});
        if (MathLib::sqrDist(q, q_projected) <= eps_pnt_out_of_plane)
            return true;
    }

    return false;
}

bool barycentricPointInTriangle(MathLib::Point3d const& p,
                                MathLib::Point3d const& a,
                                MathLib::Point3d const& b,
                                MathLib::Point3d const& c,
                                double eps_pnt_out_of_plane,
                                double eps_pnt_out_of_tri)
{
    if (std::abs(MathLib::orientation3d(p, a, b, c)) > eps_pnt_out_of_plane)
        return false;

    MathLib::Vector3 const pa(p, a);
    MathLib::Vector3 const pb(p, b);
    MathLib::Vector3 const pc(p, c);
    double const area_x_2(calcTriangleArea(a, b, c) * 2);

    double const alpha((MathLib::crossProduct(pb, pc).getLength()) / area_x_2);
    if (alpha < -eps_pnt_out_of_tri || alpha > 1 + eps_pnt_out_of_tri)
        return false;
    double const beta((MathLib::crossProduct(pc, pa).getLength()) / area_x_2);
    if (beta < -eps_pnt_out_of_tri || beta > 1 + eps_pnt_out_of_tri)
        return false;
    double const gamma(1 - alpha - beta);
    return !(gamma < -eps_pnt_out_of_tri || gamma > 1 + eps_pnt_out_of_tri);
}

bool isPointInTriangleXY(MathLib::Point3d const& p,
                         MathLib::Point3d const& a,
                         MathLib::Point3d const& b,
                         MathLib::Point3d const& c)
{
    // criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
    Eigen::Matrix2d mat;
    mat(0, 0) = b[0] - a[0];
    mat(0, 1) = c[0] - a[0];
    mat(1, 0) = b[1] - a[1];
    mat(1, 1) = c[1] - a[1];
    Eigen::Vector2d y;
    y << p[0] - a[0], p[1] - a[1];

    y = mat.partialPivLu().solve(y);

    // check if u0 and u1 fulfills the condition
    return 0 <= y[0] && y[0] <= 1 && 0 <= y[1] && y[1] <= 1 && y[0] + y[1] <= 1;
}

bool dividedByPlane(const MathLib::Point3d& a, const MathLib::Point3d& b,
                    const MathLib::Point3d& c, const MathLib::Point3d& d)
{
    for (unsigned x = 0; x < 3; ++x)
    {
        const unsigned y = (x + 1) % 3;
        const double abc =
            (b[x] - a[x]) * (c[y] - a[y]) - (b[y] - a[y]) * (c[x] - a[x]);
        const double abd =
            (b[x] - a[x]) * (d[y] - a[y]) - (b[y] - a[y]) * (d[x] - a[x]);

        if ((abc > 0 && abd < 0) || (abc < 0 && abd > 0))
            return true;
    }
    return false;
}

bool isCoplanar(const MathLib::Point3d& a, const MathLib::Point3d& b,
                const MathLib::Point3d& c, const MathLib::Point3d& d)
{
    const MathLib::Vector3 ab(a, b);
    const MathLib::Vector3 ac(a, c);
    const MathLib::Vector3 ad(a, d);

    if (ab.getSqrLength() < pow(std::numeric_limits<double>::epsilon(), 2) ||
        ac.getSqrLength() < pow(std::numeric_limits<double>::epsilon(), 2) ||
        ad.getSqrLength() < pow(std::numeric_limits<double>::epsilon(), 2))
    {
        return true;
    }

    // In exact arithmetic <ac*ad^T, ab> should be zero
    // if all four points are coplanar.
    const double sqr_scalar_triple(
        pow(MathLib::scalarProduct(MathLib::crossProduct(ac, ad), ab), 2));
    // Due to evaluating the above numerically some cancellation or rounding
    // can occur. For this reason a normalisation factor is introduced.
    const double normalisation_factor =
        (ab.getSqrLength() * ac.getSqrLength() * ad.getSqrLength());

    // tolerance 1e-11 is choosen such that
    // a = (0,0,0), b=(1,0,0), c=(0,1,0) and d=(1,1,1e-6) are considered as
    // coplanar
    // a = (0,0,0), b=(1,0,0), c=(0,1,0) and d=(1,1,1e-5) are considered as not
    // coplanar
    return (sqr_scalar_triple / normalisation_factor < 1e-11);
}

}  // end namespace MathLib
