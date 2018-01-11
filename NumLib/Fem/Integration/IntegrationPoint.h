/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/TemplateWeightedPoint.h"

namespace NumLib
{
/// Integration rule for point elements.
///
/// The integration order is not stored or used for point integration.
/// It is only needed to satisfy the common integration rule concepts.
class IntegrationPoint
{
    using WeightedPoint = MathLib::TemplateWeightedPoint<double, double, 1>;

public:
    /// IntegrationPoint constructor for given order.
    explicit IntegrationPoint(unsigned /* order */)
    {
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned /* order */)
    {
    }

    /// Return current integration order.
    unsigned getIntegrationOrder() const
    {
        return 0;
    }

    /// Return the number of sampling points.
    unsigned getNumberOfPoints() const
    {
        return 1;
    }

    /// \copydoc IntegrationGaussRegular::getWeightedPoint(unsigned) const
    WeightedPoint getWeightedPoint(unsigned igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// \copydoc IntegrationGaussRegular::getWeightedPoint(unsigned, unsigned)
    static WeightedPoint getWeightedPoint(unsigned order, unsigned igp)
    {
        (void)order;
        (void)igp;
        return WeightedPoint({{1}}, 1);
    }

    template <typename Method>
    static WeightedPoint getWeightedPoint(unsigned /*igp*/)
    {
        return WeightedPoint({{1}}, 1);
    }

    /// Get the number of integration points.
    ///
    /// \param order    the number of integration points
    /// \return the number of points.
    static unsigned getNumberOfPoints(unsigned order)
    {
        (void)order;
        return 1;
    }
};
}
