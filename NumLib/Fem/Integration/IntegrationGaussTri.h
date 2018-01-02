/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/Integration/GaussLegendreTri.h"

namespace NumLib
{

/**
 * \brief Gauss quadrature rule for triangles
 *
 * Gauss quadrature rule for triangles is originally given as
 * \f[
 *    \int F(x,y) dx dy = \int F(x(r, s), y(r, s)) j(r,s) dr ds \approx \frac{1}{2} \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w_i )
 * \f]
 *
 * To make it consistent with other elements, we rewrite the above formula as
 * \f[
 *    \int F(x,y) dx dy \approx \sum_i ( F(x(r_i, s_i), y(r_i, s_i)) w'_i )
 * \f]
 * by defining the new weight \f$ w'=\frac{1}{2} w \f$.
 */
class IntegrationGaussTri
{
    using WeightedPoint = MathLib::TemplateWeightedPoint<double, double, 2>;

public:
    /**
     * Construct this object with the given integration order
     *
     * @param order     integration order (default 2)
     */
    explicit IntegrationGaussTri(unsigned order = 2)
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
            case 1: return getWeightedPoint<MathLib::GaussLegendreTri<1> >(igp);
            case 2: return getWeightedPoint<MathLib::GaussLegendreTri<2> >(igp);
            case 3: return getWeightedPoint<MathLib::GaussLegendreTri<3> >(igp);
        }
        return WeightedPoint(std::array<double, 2>(), 0);
    }

    template <typename Method>
    static WeightedPoint
    getWeightedPoint(unsigned igp)
    {
        return WeightedPoint(Method::X[igp], 0.5 * Method::W[igp]);
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
            case 1: return MathLib::GaussLegendreTri<1>::NPoints;
            case 2: return MathLib::GaussLegendreTri<2>::NPoints;
            case 3: return MathLib::GaussLegendreTri<3>::NPoints;
        }
        return 0;
    }

private:
    unsigned _order;
    unsigned _n_sampl_pt;
};

}
