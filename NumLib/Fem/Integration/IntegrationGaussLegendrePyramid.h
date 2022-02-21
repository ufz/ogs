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

#include "BaseLib/Error.h"
#include "MathLib/Integration/GaussLegendrePyramid.h"
#include "MathLib/WeightedPoint.h"

namespace NumLib
{
/**
 * \brief Gauss-Legendre quadrature rule for pyramid
 */
class IntegrationGaussLegendrePyramid
{
public:
    /**
     * Construct this object with the given integration order
     *
     * @param order     integration order (default 2)
     */
    explicit IntegrationGaussLegendrePyramid(unsigned const order = 2)
        : _order(order)
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
    /**
     * get coordinates of a integration point
     *
     * @param igp      The integration point index
     * @return a weighted point
     */
    MathLib::WeightedPoint getWeightedPoint(unsigned const igp) const
    {
        return getWeightedPoint(getIntegrationOrder(), igp);
    }

    /**
     * get coordinates of a integration point
     *
     * @param order    the number of integration points
     * @param igp      the sampling point id
     * @return weight
     */
    static MathLib::WeightedPoint getWeightedPoint(unsigned const order,
                                                   unsigned const igp)
    {
        // For the case of order = 4, it
        // causes `assertion `rank == num_nodes' failed`
        // in the SVD decomposition in the least square extrapolation.
        // Therefore for the case of order = 4, the rule of order 3 is used.
        switch (order)
        {
            case 1:
                return getWeightedPoint<MathLib::GaussLegendrePyramid<1>>(igp);
            case 2:
                return getWeightedPoint<MathLib::GaussLegendrePyramid<2>>(igp);
            case 3:
                return getWeightedPoint<MathLib::GaussLegendrePyramid<3>>(igp);
            case 4:
                return getWeightedPoint<MathLib::GaussLegendrePyramid<3>>(igp);
        }
        OGS_FATAL("Integration order {:d} not implemented for pyramids.",
                  order);
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
        // For the case of order = 4, it
        // causes `assertion `rank == num_nodes' failed`
        // in the SVD decomposition in the least square extrapolation.
        // Therefore for the case of order = 4, the rule of order 3 is used.
        switch (order)
        {
            case 1:
                return MathLib::GaussLegendrePyramid<1>::NPoints;
            case 2:
                return MathLib::GaussLegendrePyramid<2>::NPoints;
            case 3:
                return MathLib::GaussLegendrePyramid<3>::NPoints;
            case 4:
                return MathLib::GaussLegendrePyramid<3>::NPoints;
        }
        OGS_FATAL("Integration order {:d} not implemented for pyramids.",
                  order);
    }

private:
    unsigned _order;
    unsigned _n_sampl_pt{0};
};

}  // namespace NumLib
