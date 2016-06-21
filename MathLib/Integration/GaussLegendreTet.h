/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GAUSSLEGENDRETET_H_
#define GAUSSLEGENDRETET_H_

#include <array>
#ifdef _MSC_VER
#include "mathlib_export.h"
#else
#define MATHLIB_EXPORT
#endif

namespace MathLib
{

/// Gauss-Legendre quadrature on tetrahedrals
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendreTet {
    static MATHLIB_EXPORT const unsigned Order = ORDER;
    static MATHLIB_EXPORT const unsigned NPoints = ORDER;
    static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendreTet<2> {
    static MATHLIB_EXPORT const unsigned Order = 2;
    static MATHLIB_EXPORT const unsigned NPoints = 5;
    static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendreTet<3> {
    static MATHLIB_EXPORT const unsigned Order = 3;
    static MATHLIB_EXPORT const unsigned NPoints = 15;
    static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

}

#endif //GAUSSLEGENDRETET_H_
