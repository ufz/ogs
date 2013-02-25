/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-04
 * \brief  Definition of the Color class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef COLOR_H_
#define COLOR_H_

#include "TemplatePoint.h"

#include <cstdlib>
#include <fstream>
#include <list>
#include <map>

namespace GeoLib
{
typedef MathLib::TemplatePoint<unsigned char> Color;

/// Returns a random RGB colour.
Color* getRandomColor();

/// Reads a color-lookup-table from a file and writes it to a map.
int readColorLookupTable(std::map<std::string, GeoLib::Color*> &colors, const std::string &filename);

/// Uses a color-lookup-table (in form of a map) to return a colour for a specified name. If the name is not
/// in the colortable a new entry is created with the new name and a random colour.
const Color* getColor(const std::string &id, std::map<std::string, GeoLib::Color*> &colors);

/// Convenience function to use the getColor method with numbers as identifiers.
const Color* getColor(double id, std::map<std::string, GeoLib::Color*> &colors);
} // namespace

#endif /* COLOR_H_ */
