/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on Oct 29, 2021, 12:08 AM
 *
 */
#include "IntegrationGaussLegendrePrism.h"

#include <cassert>

namespace NumLib
{
template <int OrderGaussLegendreTri, int OrderGaussLegendre>
MathLib::WeightedPoint getWeightedPointConcrete(unsigned const igp)
{
    assert(igp < MathLib::GaussLegendreTri<OrderGaussLegendreTri>::NPoints *
                     OrderGaussLegendre);

    const unsigned gp_r =
        igp % MathLib::GaussLegendreTri<OrderGaussLegendreTri>::NPoints;
    const unsigned gp_t =
        (unsigned)(igp /
                   MathLib::GaussLegendreTri<OrderGaussLegendreTri>::NPoints);
    std::array<double, 3> rst;
    rst[0] = MathLib::GaussLegendreTri<OrderGaussLegendreTri>::X[gp_r][0];
    rst[1] = MathLib::GaussLegendreTri<OrderGaussLegendreTri>::X[gp_r][1];
    rst[2] = MathLib::GaussLegendre<OrderGaussLegendre>::X[gp_t];
    double const w = MathLib::GaussLegendreTri<OrderGaussLegendreTri>::W[gp_r] *
                     0.5 * MathLib::GaussLegendre<OrderGaussLegendre>::W[gp_t];

    return MathLib::WeightedPoint(rst, w);
}

void IntegrationGaussLegendrePrism::setIntegrationOrder(unsigned const order)
{
    _order = order;
    _n_sampl_pt = getNumberOfPoints(_order);
}

MathLib::WeightedPoint IntegrationGaussLegendrePrism::getWeightedPoint(
    unsigned const order, unsigned const igp)
{
    return (order < 2) ? getWeightedPointConcrete<2, 2>(igp)
                       : getWeightedPointConcrete<4, 3>(igp);
}

unsigned IntegrationGaussLegendrePrism::getNumberOfPoints(unsigned const order)
{
    return (order < 2) ? MathLib::GaussLegendreTri<2>::NPoints * 2
                       : MathLib::GaussLegendreTri<4>::NPoints * 3;
}
}  // namespace NumLib
