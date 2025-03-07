/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

#include "mathlib_export.h"

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
    static MATHLIB_EXPORT const unsigned NPoints = 14;
    static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendreTet<4>
{
    static MATHLIB_EXPORT const unsigned Order = 4;
    static MATHLIB_EXPORT const unsigned NPoints = 20;
    static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
    static MATHLIB_EXPORT const double W[NPoints];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler but is necessary.
template <>
const std::array<std::array<double, 3>, GaussLegendreTet<1>::NPoints>
    GaussLegendreTet<1>::X;
template <>
double const GaussLegendreTet<1>::W[1];
#endif
}  // namespace MathLib
