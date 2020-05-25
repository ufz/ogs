/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ColorLookupTable.h"

#include <limits>

namespace DataHolderLib
{
ColorLookupTable::ColorLookupTable()
    : range_(
          std::make_pair<double, double>(std::numeric_limits<double>::lowest(),
                                         std::numeric_limits<double>::max()))

{}

void ColorLookupTable::setTableRange(double min, double max)
{
    if (min < max)
    {
        range_ = std::make_pair(min, max);
    }
}

void ColorLookupTable::setColor(double id, DataHolderLib::Color const& color)
{
    if ((id > range_.first) && (id < range_.second))
    {
        lut_.emplace_back(id, color, "");
    }
}

void ColorLookupTable::setColor(std::string const& name, DataHolderLib::Color const& color)
{
    lut_.emplace_back(0, color, name);
}

}  // namespace DataHolderLib
