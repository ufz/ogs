/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOMETRIC_BASICS_H_
#define GEOMETRIC_BASICS_H_

#include <cstddef>

namespace MathLib
{

template <typename T, std::size_t DIM> class TemplatePoint;
typedef MathLib::TemplatePoint<double,3> Point3d;

/**
 * Calculates the volume of a tetrahedron.
 * The formula is V=1/6*|a(b x c)| with a=x1->x2, b=x1->x3 and c=x1->x4.
 */
double calcTetrahedronVolume(MathLib::Point3d const& x1,
                             MathLib::Point3d const& x2,
                             MathLib::Point3d const& x3,
                             MathLib::Point3d const& x4);

}  // end namespace MathLib

#endif
