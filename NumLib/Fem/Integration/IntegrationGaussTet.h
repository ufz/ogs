/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/Integration/GaussLegendreTet.h"

namespace NumLib
{

/**
 * \brief Gauss quadrature rule for tetrahedrals
 */
class IntegrationGaussTet
{
    using WeightedPoint = MathLib::TemplateWeightedPoint<double, double, 3>;

public:
    /**
     * Construct this object with the given integration order
     *
     * @param order     integration order (default 2)
     */
    explicit IntegrationGaussTet(unsigned order = 2)
    : _order(order), _n_sampl_pt(0)
    {
        this->setIntegrationOrder(order);
    }

    /// Change the integration order.
    void setIntegrationOrder(unsigned order)
    {
        _order = order;
        _n_sampl_pt = getNumberOfPoints(order);
    }

    /// return current integration order.
    unsigned getIntegrationOrder() const {return _order;}

    /// return the number of sampling points
    unsigned getNumberOfPoints() const {return _n_sampl_pt;}

    /**
     * get coordinates of a integration point
     *
     * @param igp      The integration point index
     * @return a weighted point
     */
    WeightedPoint getWeightedPoint(unsigned igp) const
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
    static WeightedPoint
    getWeightedPoint(unsigned order, unsigned igp)
    {
        switch (order)
        {
            case 1: return getWeightedPoint<MathLib::GaussLegendreTet<1> >(igp);
            case 2: return getWeightedPoint<MathLib::GaussLegendreTet<2> >(igp);
            case 3: return getWeightedPoint<MathLib::GaussLegendreTet<3> >(igp);
        }
        return WeightedPoint(std::array<double, 3>(), 0);
    }

    template <typename Method>
    static WeightedPoint
    getWeightedPoint(unsigned igp)
    {
        return WeightedPoint(Method::X[igp], Method::W[igp]);
    }


    /**
     * get the number of integration points
     *
     * @param order    the number of integration points
     * @return the number of points
     */
    static unsigned
    getNumberOfPoints(unsigned order)
    {
        switch (order)
        {
            case 1: return MathLib::GaussLegendreTet<1>::NPoints;
            case 2: return MathLib::GaussLegendreTet<2>::NPoints;
            case 3: return MathLib::GaussLegendreTet<3>::NPoints;
        }
        return 0;
    }

private:
    unsigned _order;
    unsigned _n_sampl_pt;
};

}
