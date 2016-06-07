/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeometricBasics.h"
#include "Point3d.h"
#include "Vector3.h"

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
        if (!(d0_sign == (d4>=0) || std::abs(d4) < eps))
            return false;
        return true;
    }
    return false;
}

}  // end namespace MathLib
