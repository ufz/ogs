/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
struct GaussLegendre { };

template <>
struct GaussLegendre<1>
{
    static const unsigned Order;
    static const double X[1];
    static const double W[1];
};
const unsigned GaussLegendre<1>::Order = 1;
const double GaussLegendre<1>::X[1] = {0.};
const double GaussLegendre<1>::W[1] = {2.};

template<>
struct GaussLegendre<2>
{
    static const unsigned Order;
    static const double X[2];
    static const double W[2];
};
unsigned const GaussLegendre<2>::Order = 2;
double const GaussLegendre<2>::X[2] = {0.577350269189626, -0.577350269189626};
double const GaussLegendre<2>::W[2] = {1., 1.};

template<>
struct GaussLegendre<3>
{
    static const unsigned Order;
    static const double X[3];
    static const double W[3];
};
unsigned const GaussLegendre<3>::Order = 3;
double const GaussLegendre<3>::X[3] = {0.774596669241483, 0., -0.774596669241483};
double const GaussLegendre<3>::W[3] = {5./9, 8./9, 5./9};

template<>
struct GaussLegendre<4>
{
    static const unsigned Order;
    static const double X[4];
    static const double W[4];
};
unsigned const GaussLegendre<4>::Order = 4;
double const GaussLegendre<4>::X[4] =
        {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
double const GaussLegendre<4>::W[4] =
        { 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454};

}   // namespace MathLib

#endif // GAUSSLEGENDRE_H_

