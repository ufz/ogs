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

namespace
{
template <int OrderGaussLegendreTri, int OrderGaussLegendre>
MathLib::WeightedPoint getWeightedPointConcrete(unsigned const igp)
{
    using GL = MathLib::GaussLegendre<OrderGaussLegendre>;
    using GLT = MathLib::GaussLegendreTri<OrderGaussLegendreTri>;

    assert(igp < GLT::NPoints * OrderGaussLegendre);

    const unsigned gp_r = igp % GLT::NPoints;
    const unsigned gp_t = igp / GLT::NPoints;

    std::array<double, 3> rst{GLT::X[gp_r][0], GLT::X[gp_r][1], GL::X[gp_t]};

    double const weight = GLT::W[gp_r] * 0.5 * GL::W[gp_t];

    return MathLib::WeightedPoint(rst, weight);
}
}  // namespace

namespace NumLib
{

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
