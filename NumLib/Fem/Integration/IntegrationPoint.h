/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_INTEGRATIONPOINT_H_
#define NUMLIB_INTEGRATIONPOINT_H_

namespace NumLib
{
/// Integration rule for point elements.
///
/// The integration order is not stored or used for point integration.
/// It is only needed to satisfy the common integration rule concepts.
class IntegrationPoint
{
    typedef MathLib::TemplateWeightedPoint<double, double, 1> WeightedPoint;

public:
    /// IntegrationPoint constructor for given order.
    explicit IntegrationPoint(std::size_t /* order */)
    {
    }

    /// Change the integration order.
    void setIntegrationOrder(std::size_t /* order */)
    {
    }

    /// Return current integration order.
    std::size_t getIntegrationOrder() const
    {
        return 0;
    }

    /// Return the number of sampling points.
    std::size_t getNumberOfPoints() const
    {
        return 1;
    }

    /// Get coordinates of a integration point.
    ///
    /// \param igp      The integration point index.
    /// \return a weighted point.
    WeightedPoint getWeightedPoint(std::size_t igp)
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /// Get coordinates of a integration point.
    ///
    /// \param order    the number of integration points.
    /// \param igp      the sampling point id.
    /// \return weight
    static WeightedPoint getWeightedPoint(std::size_t /* order */,
                                          std::size_t /*igp*/)
    {
        return WeightedPoint({{1}}, 1);
    }

    template <typename Method>
    static WeightedPoint getWeightedPoint(std::size_t /*igp*/)
    {
        return WeightedPoint({{1}}, 1);
    }

    /// Get the number of integration points.
    ///
    /// \param order    the number of integration points
    /// \return the number of points.
    static std::size_t getNumberOfPoints(std::size_t /* order */)
    {
        return 1;
    }
};
}

#endif  // NUMLIB_INTEGRATIONPOINT_H_
