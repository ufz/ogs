/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef INTEGRATIONGAUSSREGULAR_H_
#define INTEGRATIONGAUSSREGULAR_H_

#include <cmath>
#include <array>

#include "MathLib/Integration/GaussLegendre.h"
#include "MathLib/TemplateWeightedPoint.h"


namespace NumLib
{

/// Gauss quadrature rule for regular shape elements: line, quad and hex.
///
/// \tparam N_DIM    Spatial dimension
template <std::size_t N_DIM>
class IntegrationGaussRegular
{
    typedef typename MathLib::TemplateWeightedPoint<double, double, N_DIM>
        WeightedPoint;
public:
    /// Create IntegrationGaussRegular of the given Gauss-Legendre integration
    /// order.
    ///
    /// @param order     integration order (default 2)
    explicit IntegrationGaussRegular(std::size_t order = 2)
    : _order(order), _n_sampl_pt(0)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(std::size_t order)
    {
        this->_n_sampl_pt = static_cast<std::size_t>(std::pow(order, N_DIM));
        this->_order = order;
    }

    /// return current integration order.
    std::size_t getIntegrationOrder() const {return _order;}

    /// return the number of sampling points
    std::size_t getNPoints() const {return _n_sampl_pt;}

    /// Get coordinates of the integration point.
    ///
    /// @param igp       The integration point index
    /// @return a weighted point
    WeightedPoint
    getWeightedPoint(std::size_t igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// Get position indexes in r-s-t axis.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return  a tuple of position indexes
    static std::array<std::size_t, N_DIM> getPositionIndices(std::size_t order, std::size_t igp);

    /// Get coordinates of the integration point.
    ///
    /// @param order     The number of integration points
    /// @param igp       The integration point index
    /// @return a weighted point
    static WeightedPoint
    getWeightedPoint(std::size_t order, std::size_t igp);

private:
    /// Computes weighted point using given integration method.
    ///
    /// \tparam Method  Integration method to use.
    /// \param  pos     Point indices computed by getPositionIndices.
    template <typename Method>
    static
    WeightedPoint
    getWeightedPoint(std::array<std::size_t, N_DIM> const& pos);

private:
    std::size_t _order;
    std::size_t _n_sampl_pt;
};

} // NumLib

#include "IntegrationGaussRegular-impl.h"

#endif //INTEGRATIONGAUSSREGULAR_H_
