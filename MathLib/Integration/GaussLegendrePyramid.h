/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

namespace MathLib
{

/// Gauss-Legendre quadrature on pyramid
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendrePyramid {
    static const unsigned Order = ORDER;
    static const unsigned NPoints = ORDER;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

template <>
struct GaussLegendrePyramid<2> {
    static const unsigned Order = 2;
    static const unsigned NPoints = 5;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

template <>
struct GaussLegendrePyramid<3> {
    static const unsigned Order = 3;
    static const unsigned NPoints = 13;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler but is necessary.
template <>
const std::array<std::array<double, 3>, GaussLegendrePyramid<1>::NPoints>
    GaussLegendrePyramid<1>::X;
template <>
double const GaussLegendrePyramid<1>::W[1];
#endif
}
