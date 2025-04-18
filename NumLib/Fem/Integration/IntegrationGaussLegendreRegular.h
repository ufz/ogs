/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <cmath>

#include "BaseLib/Error.h"
#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/WeightedPoint.h"

namespace NumLib
{
/// Gauss-Legendre quadrature rule for regular shape elements: line, quad and
/// hex.
///
/// \tparam N_DIM    Spatial dimension
template <unsigned N_DIM>
class IntegrationGaussLegendreRegular
{
public:
    /// Create IntegrationGaussLegendreRegular of the given Gauss-Legendre
    /// integration order.
    ///
    /// @param order     integration order (default 2)
    explicit IntegrationGaussLegendreRegular(unsigned order = 2) : _order(order)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned order)
    {
        this->_n_sampl_pt = static_cast<unsigned>(std::pow(order, N_DIM));
        this->_order = order;
    }

    /// return current integration order.
    unsigned getIntegrationOrder() const { return _order; }

    /// return the number of sampling points
    unsigned getNumberOfPoints() const { return _n_sampl_pt; }

    /// Get coordinates of the integration point.
    ///
    /// @param igp       The integration point index
    /// @return a weighted point
    MathLib::WeightedPoint getWeightedPoint(unsigned const igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// Get position indexes in r-s-t axis.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return  a tuple of position indexes
    static std::array<unsigned, N_DIM> getPositionIndices(unsigned order,
                                                          unsigned igp);

    /// Get coordinates of the integration point.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return a weighted point
    static MathLib::WeightedPoint getWeightedPoint(unsigned const order,
                                                   unsigned const igp);

private:
    /// Computes weighted point using given integration method.
    ///
    /// \tparam Method  Integration method to use.
    /// \param  pos     Point indices computed by getPositionIndices.
    template <typename Method>
    static MathLib::WeightedPoint getWeightedPoint(
        std::array<unsigned, N_DIM> const& pos);

private:
    unsigned _order;
    unsigned _n_sampl_pt{0};
};

}  // namespace NumLib

#include "IntegrationGaussLegendreRegular-impl.h"
