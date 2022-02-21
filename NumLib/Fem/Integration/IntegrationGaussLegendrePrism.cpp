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

#include "BaseLib/Error.h"

namespace
{
template <int OrderGaussLegendreTri, int OrderGaussLegendre>
constexpr unsigned getNumberOfPointsConcrete()
{
    return MathLib::GaussLegendreTri<OrderGaussLegendreTri>::NPoints *
           OrderGaussLegendre;
}

template <int OrderGaussLegendreTri, int OrderGaussLegendre>
MathLib::WeightedPoint getWeightedPointConcrete(unsigned const igp)
{
    using GL = MathLib::GaussLegendre<OrderGaussLegendre>;
    using GLT = MathLib::GaussLegendreTri<OrderGaussLegendreTri>;

    assert(igp < (getNumberOfPointsConcrete<OrderGaussLegendreTri,
                                            OrderGaussLegendre>()));

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
    // Note: These cases must correspod strictly to the logic in
    // getNumberOfPoints()!
    switch (order)
    {
        case 1:
            return getWeightedPointConcrete<1, 1>(igp);
        case 2:
            return getWeightedPointConcrete<2, 2>(igp);
        case 3:
            // The combination <4, 3> has been chosen to allow extrapolation
            // from the set of integration points on Prism13 elements for the
            // third integration order. <3, 3> would not be sufficient.
            return getWeightedPointConcrete<4, 3>(igp);
        case 4:
            return getWeightedPointConcrete<4, 4>(igp);
        default:
            OGS_FATAL(
                "Integration order {} not supported for integration on prisms.",
                order);
    }
}

unsigned IntegrationGaussLegendrePrism::getNumberOfPoints(unsigned const order)
{
    // Note: These cases must correspod strictly to the logic in
    // getWeightedPoint()!
    switch (order)
    {
        case 1:
            return getNumberOfPointsConcrete<1, 1>();
        case 2:
            return getNumberOfPointsConcrete<2, 2>();
        case 3:
            // The combination <4, 3> has been chosen to allow extrapolation
            // from the set of integration points on Prism13 elements for the
            // third integration order. <3, 3> would not be sufficient.
            return getNumberOfPointsConcrete<4, 3>();
        case 4:
            return getNumberOfPointsConcrete<4, 4>();
        default:
            OGS_FATAL(
                "Integration order {} not supported for integration on prisms.",
                order);
    }
}
}  // namespace NumLib
