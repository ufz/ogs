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

}  // end namespace MathLib
