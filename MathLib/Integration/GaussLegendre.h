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


#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include "WeightedSum.h"

namespace MathLib
{

/// Gauss-Legendre quadrature method
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendre {
    static const unsigned Order = ORDER;
    static const double X[Order];
    static const double W[Order];
};

#ifndef _MSC_VER  // The following explicit instantatiation declaration does not
                  // compile on that particular compiler.
template <>
double const GaussLegendre<1>::X[1];
template <>
double const GaussLegendre<1>::W[1];
template <>
double const GaussLegendre<2>::X[2];
template <>
double const GaussLegendre<2>::W[2];
template <>
double const GaussLegendre<3>::X[3];
template <>
double const GaussLegendre<3>::W[3];
template <>
double const GaussLegendre<4>::X[4];
template <>
double const GaussLegendre<4>::W[4];
#endif

}  // namespace MathLib

#endif // GAUSSLEGENDRE_H_

