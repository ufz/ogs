/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "WeightedPoint.h"

#include <ostream>

namespace MathLib
{
// Used to output helpful information if the unit tests fail.
std::ostream& operator<<(std::ostream& os, MathLib::WeightedPoint const& wp)
{
    auto const dim = wp.getDimension();
    os << "WP[" << dim << "D]{{";
    for (std::size_t comp = 0; comp < 3; ++comp)
    {
        if (comp != 0)
        {
            os << ", ";
        }
        os << wp[comp];
    }
    os << "}, weight=" << wp.getWeight() << '}';
    return os;
}
}  // namespace MathLib
