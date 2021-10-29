/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/Integration/GaussLegendreTri.h"
#include "MathLib/WeightedPoint.h"

namespace NumLib
{
/**
 * \brief Gauss-Legendre quadrature rule for prisms
 */
class IntegrationGaussLegendrePrism
{
public:
    /**
     * Construct this object with the given integration order
     *
     * @param order     integration order (default 2)
     */
    explicit IntegrationGaussLegendrePrism(unsigned const order)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned const order);

    /// return current integration order.
    unsigned getIntegrationOrder() const { return _order; }

    /// return the number of sampling points
    unsigned getNumberOfPoints() const { return _n_sampl_pt; }

    // clang-format off
    /// \copydoc NumLib::IntegrationGaussLegendreRegular::getWeightedPoint(unsigned) const
    // clang-format on
    MathLib::WeightedPoint getWeightedPoint(unsigned const igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    // clang-format off
    /// \copydoc NumLib::IntegrationGaussLegendreRegular::getWeightedPoint(unsigned, unsigned)
    // clang-format on
    static MathLib::WeightedPoint getWeightedPoint(unsigned const order,
                                                   unsigned const igp);

    template <typename Method>
    static MathLib::WeightedPoint getWeightedPoint(unsigned const igp)
    {
        return MathLib::WeightedPoint(Method::X[igp], Method::W[igp]);
    }

    /**
     * get the number of integration points
     *
     * @param order    the number of integration points
     * @return the number of points
     */
    static unsigned getNumberOfPoints(unsigned const order);

private:
    unsigned _order{2};
    unsigned _n_sampl_pt{0};
};

}  // namespace NumLib
