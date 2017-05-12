/**
 * \file
 * \date   2015-01-16
 * \brief  Definition of the Point3d class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>

#include "mathlib_export.h"

#include "TemplatePoint.h"
#include "MathTools.h"

namespace MathLib
{
using Point3d = MathLib::TemplatePoint<double, 3>;

extern MATHLIB_EXPORT const Point3d ORIGIN;
/**
 * rotation of points
 * @param mat a rotation matrix
 * @param p   a point to be transformed
 * @return a rotated point
 */
template <typename MATRIX>
inline MathLib::Point3d operator*(MATRIX const& mat, MathLib::Point3d const& p)
{
    MathLib::Point3d new_p;
    for (std::size_t i(0); i<3; ++i) {
        for (std::size_t j(0); j<3; ++j) {
            new_p[i] += mat(i,j)*p[j];
        }
    }
    return new_p;
}

/** Computes the squared dist between the two points p0 and p1.
 */
inline
double sqrDist(MathLib::Point3d const& p0, MathLib::Point3d const& p1)
{
    const double v[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    return MathLib::scalarProduct<double,3>(v,v);
}

/// Computes the squared distance between the orthogonal projection of the two
/// points \c p0 and \c p1 onto the \f$xy\f$-plane.
inline
double sqrDist2d(MathLib::Point3d const& p0, MathLib::Point3d const& p1)
{
    return (p0[0]-p1[0])*(p0[0]-p1[0]) + (p0[1]-p1[1])*(p0[1]-p1[1]);
}

} // end namespace MathLib
