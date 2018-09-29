/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <tuple>
#include <utility>
#include <vector>

#include "Applications/DataHolderLib/Color.h"

namespace DataHolderLib
{

/// Interpolation methods
enum class LUTType {
    NONE = 0,
    LINEAR = 1,
    EXPONENTIAL = 2,
    SIGMOID = 3 // not yet implemented
};

/**
 * A colour look-up table stored as a vector containing for each entry
 *    id     - value for which the colour is stored
 *    colour - RGBA value
 *    name   - a name referencing the colour (such as a stratigraphic layer)
 * Colours from the table can then be accessed using either id or name.
 */
class ColorLookupTable
{
public:
    ColorLookupTable();

    void setColor(double id, DataHolderLib::Color const& color);

    void setColor(std::string const& name, DataHolderLib::Color const& color);

    DataHolderLib::LUTType getInterpolationType() const { return _type; }

    void setInterpolationType(LUTType type) { _type = type; }

    std::size_t size() const { return _lut.size(); }

    std::pair<double, double> getTableRange() const { return _range; }

    void setTableRange(double min, double max);

    std::tuple<double, Color, std::string> const& operator[](std::size_t i) const
    {
        assert (i < _lut.size());
        return _lut[i];
    }

private:
    std::vector< std::tuple<double, Color, std::string> > _lut;
    LUTType _type;
    std::pair<double, double> _range;
};

} // namespace
