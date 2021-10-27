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
    explicit IntegrationGaussLegendrePrism(unsigned order = 2)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned const order)
    {
        _order = order;
        _n_sampl_pt = getNumberOfPoints(_order);
    }

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
                                                   unsigned const igp)
    {
        if (order < 3)
        {
            const unsigned gp_r = igp % 3;
            const auto gp_t = (unsigned)(igp / 3);
            std::array<double, 3> rst;
            rst[0] = MathLib::GaussLegendreTri<2>::X[gp_r][0];
            rst[1] = MathLib::GaussLegendreTri<2>::X[gp_r][1];
            rst[2] = MathLib::GaussLegendre<2>::X[gp_t];
            double w = MathLib::GaussLegendreTri<2>::W[gp_r] * 0.5 *
                       MathLib::GaussLegendre<2>::W[gp_t];
            return MathLib::WeightedPoint(rst, w);
        }
        const unsigned gp_r = igp % MathLib::GaussLegendreTri<4>::NPoints;
        const auto gp_t =
            (unsigned)(igp / MathLib::GaussLegendreTri<4>::NPoints);
        std::array<double, 3> rst;
        rst[0] = MathLib::GaussLegendreTri<4>::X[gp_r][0];
        rst[1] = MathLib::GaussLegendreTri<4>::X[gp_r][1];
        rst[2] = MathLib::GaussLegendre<3>::X[gp_t];
        double w = MathLib::GaussLegendreTri<4>::W[gp_r] * 0.5 *
                   MathLib::GaussLegendre<3>::W[gp_t];
        return MathLib::WeightedPoint(rst, w);
    }

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
    static unsigned getNumberOfPoints(unsigned const order)
    {
        return (order < 3) ? 6 : 21;
    }

private:
    unsigned _order{2};
    unsigned _n_sampl_pt{0};
};

}  // namespace NumLib
