// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SpatialPosition.h"

#include <iomanip>
#include <limits>
#include <ostream>

namespace ParameterLib
{
std::ostream& operator<<(std::ostream& os, SpatialPosition const& pos)
{
    auto const default_precision = os.precision();
    os << "{ node: ";
    if (auto const node_id = pos.getNodeID())
    {
        os << *node_id;
    }
    else
    {
        os << "not set";
    }

    os << ", element: ";
    if (auto const element_id = pos.getElementID())
    {
        os << *element_id;
    }
    else
    {
        os << "not set";
    }

    os << ", coordinates: ";
    if (auto const coordinates = pos.getCoordinates())
    {
        os << std::setprecision(std::numeric_limits<double>::digits10 + 1)
           << "(" << *coordinates << ")"
           << std::setprecision(default_precision);
    }
    else
    {
        os << "not set";
    }
    os << " }";
    return os;
}
}  // namespace ParameterLib
