// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
