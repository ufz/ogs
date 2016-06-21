/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include "WeightedSum.h"
#ifdef _MSC_VER
#include "mathlib_export.h"
#else
#define MATHLIB_EXPORT
#endif

namespace MathLib
{

/// Gauss-Legendre quadrature method
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendre {
	static MATHLIB_EXPORT const unsigned Order = ORDER;
	static MATHLIB_EXPORT const double X[Order];
	static MATHLIB_EXPORT const double W[Order];
};

}   // namespace MathLib

#endif // GAUSSLEGENDRE_H_

