// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ColorLookupTable.h"

#include <limits>

namespace DataHolderLib
{
ColorLookupTable::ColorLookupTable()
    : _range(
          std::make_pair<double, double>(std::numeric_limits<double>::lowest(),
                                         std::numeric_limits<double>::max()))

{
}

void ColorLookupTable::setTableRange(double min, double max)
{
    if (min < max)
    {
        _range = std::make_pair(min, max);
    }
}

void ColorLookupTable::setColor(double id, DataHolderLib::Color const& color)
{
    if ((id > _range.first) && (id < _range.second))
    {
        _lut.emplace_back(id, color, "");
    }
}

void ColorLookupTable::setColor(std::string const& name,
                                DataHolderLib::Color const& color)
{
    _lut.emplace_back(0, color, name);
}

}  // namespace DataHolderLib
