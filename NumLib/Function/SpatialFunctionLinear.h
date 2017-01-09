/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SPATIALFUNCTIONLINEAR_H_
#define SPATIALFUNCTIONLINEAR_H_

#include "MathLib/LinearFunction.h"
#include "TemplateSpatialFunction.h"

namespace NumLib
{

/// Representation of linear functions f: R^3 -> R; f({x, y, z}) = a + bx + cy + dz
typedef TemplateSpatialFunction<MathLib::LinearFunction<double, 3> > SpatialFunctionLinear;

}

#endif /* SPATIALFUNCTIONLINEAR_H_ */
