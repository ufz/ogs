/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GAUSSLEGENDRETRI_H_
#define GAUSSLEGENDRETRI_H_

#include <array>
#ifdef _MSC_VER
#include "mathlib_export.h"
#else
#define MATHLIB_EXPORT
#endif

namespace MathLib
{

/// Gauss-Legendre quadrature on triangles
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendreTri {
    static const unsigned Order = ORDER;
    static const unsigned NPoints = ORDER;
	static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendreTri<2> {
    static const unsigned Order = 2;
    static const unsigned NPoints = 3;
	static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendreTri<3> {
    static const unsigned Order = 3;
    static const unsigned NPoints = 4;
	static MATHLIB_EXPORT const std::array<std::array<double, 2>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler but is necessary.
template <>
const std::array<std::array<double, 2>, GaussLegendreTri<1>::NPoints>
    GaussLegendreTri<1>::X;
template <>
double const GaussLegendreTri<1>::W[1];
#endif
}

#endif //GAUSSLEGENDRETRI_H_
