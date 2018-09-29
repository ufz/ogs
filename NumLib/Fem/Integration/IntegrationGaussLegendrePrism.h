/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/Integration/GaussLegendreTri.h"
#include "MathLib/TemplateWeightedPoint.h"

namespace NumLib
{
/**
 * \brief Gauss-Legendre quadrature rule for prisms
 */
class IntegrationGaussLegendrePrism
{
    using WeightedPoint = MathLib::TemplateWeightedPoint<double, double, 3>;

public:
    /**
     * Construct this object with the given integration order
     *
     * @param order     integration order (default 2)
     */
    explicit IntegrationGaussLegendrePrism(unsigned order = 2)
        : _order(2), _n_sampl_pt(0)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned /*order*/)
    {
        _order = 2;  // fixed
        _n_sampl_pt = getNumberOfPoints(_order);
    }

    /// return current integration order.
    unsigned getIntegrationOrder() const { return _order; }

    /// return the number of sampling points
    unsigned getNumberOfPoints() const { return _n_sampl_pt; }

    /// \copydoc NumLib::IntegrationGaussLegendreRegular::getWeightedPoint(unsigned) const
    WeightedPoint getWeightedPoint(unsigned igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// \copydoc NumLib::IntegrationGaussLegendreRegular::getWeightedPoint(unsigned, unsigned)
    static WeightedPoint getWeightedPoint(unsigned order, unsigned igp)
    {
        (void)order;
        const unsigned gp_r = igp % 3;
        const auto gp_t = (unsigned)(igp / 3);
        std::array<double, 3> rst;
        rst[0] = MathLib::GaussLegendreTri<2>::X[gp_r][0];
        rst[1] = MathLib::GaussLegendreTri<2>::X[gp_r][1];
        rst[2] = MathLib::GaussLegendre<2>::X[gp_t];
        double w = MathLib::GaussLegendreTri<2>::W[gp_r] * 0.5 *
                   MathLib::GaussLegendre<2>::W[gp_t];
        return WeightedPoint(rst, w);
    }

    template <typename Method>
    static WeightedPoint getWeightedPoint(unsigned igp)
    {
        return WeightedPoint(Method::X[igp], Method::W[igp]);
    }

    /**
     * get the number of integration points
     *
     * @param order    the number of integration points
     * @return the number of points
     */
    static unsigned getNumberOfPoints(unsigned order)
    {
        if (order == 2)
            return 6;
        return 0;
    }

private:
    unsigned _order;
    unsigned _n_sampl_pt;
};

}  // namespace NumLib
