/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-04
 * \brief  Definition of the Color class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#pragma once

#include <array>
#include <map>
#include <string>

namespace DataHolderLib
{
using Color = std::array<unsigned char, 4>;

Color createColor(unsigned char r, unsigned char g, unsigned char b);

Color createColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a);

/// Returns a random RGB colour.
Color getRandomColor();

/// Uses a color-lookup-table (in form of a map) to return a colour for a specified name. If the name is not
/// in the colortable a new entry is created with the new name and a random colour.
Color const getColor(const std::string &id, std::map<std::string, DataHolderLib::Color> &colors);

/// Convenience function to use the getColor method with numbers as identifiers.
Color const getColor(double id, std::map<std::string, DataHolderLib::Color> &colors);
} // namespace
