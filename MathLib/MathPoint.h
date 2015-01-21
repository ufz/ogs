/**
 * \file
 * \date   2015-01-16
 * \brief  Definition of the MathPoint class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHPOINT_H_
#define MATHPOINT_H_

#include <limits>

#include "TemplatePoint.h"

namespace MathLib
{
typedef MathLib::TemplatePoint<double,3> MathPoint;
} // end namespace MathLib

/**
 * lexicographic comparison of points
 */
bool operator<= (MathLib::MathPoint const & p0, MathLib::MathPoint const & p1);

/**
 * lexicographical comparison of points taking an epsilon into account
 * @param p0 first input MathPoint
 * @param p1 first input MathPoint
 * @param tol tolerance (if in the comparison operation the property fabs(p0[k] - p1[k]) < tol
 *     holds for the k-th coordinate the points are assumed the be equal in this coordinate)
 * @return true, if p0 is lexicographically smaller than p1
 */
bool lessEq(const MathLib::MathPoint& p0,
            const MathLib::MathPoint& p1,
            double tol = std::numeric_limits<double>::epsilon());

#endif /* MATHPOINT_H_ */

