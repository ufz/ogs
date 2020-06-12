/**
 * \file
 * \date Sep 3, 2013
 * \brief Weighted point class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "TemplatePoint.h"

namespace MathLib
{

template <typename FP_T, typename W_T, std::size_t DIM>
class TemplateWeightedPoint : public TemplatePoint<FP_T, DIM>
{
public:
    TemplateWeightedPoint(std::array<FP_T, DIM> const& x, W_T weight) :
        TemplatePoint<FP_T, DIM>(x), weight_(weight)
    {}

    W_T getWeight() const
    {
        return weight_;
    }

private:
    W_T const weight_;
};

using WeightedPoint1D = TemplateWeightedPoint<double, double, 1>;
using WeightedPoint2D = TemplateWeightedPoint<double, double, 2>;
using WeightedPoint3D = TemplateWeightedPoint<double, double, 3>;

} // end namespace MathLib
