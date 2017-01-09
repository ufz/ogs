/**
 * @file TemplateWeightedPoint.h
 * @date Sep 3, 2013
 * @brief Weighted point class.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef TEMPLATEWEIGHTEDPOINT_H_
#define TEMPLATEWEIGHTEDPOINT_H_

#include "TemplatePoint.h"

namespace MathLib
{

template <typename FP_T, typename W_T, std::size_t DIM>
class TemplateWeightedPoint : public TemplatePoint<FP_T, DIM>
{
public:
    TemplateWeightedPoint(std::array<FP_T, DIM> const& x, W_T weight) :
        TemplatePoint<FP_T, DIM>(x), _weight(weight)
    {}

    W_T getWeight() const
    {
        return _weight;
    }

private:
    W_T const _weight;
};

typedef TemplateWeightedPoint<double, double, 1> WeightedPoint1D;
typedef TemplateWeightedPoint<double, double, 2> WeightedPoint2D;
typedef TemplateWeightedPoint<double, double, 3> WeightedPoint3D;

} // end namespace MathLib


#endif /* TEMPLATEWEIGHTEDPOINT_H_ */
