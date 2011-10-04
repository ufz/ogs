/**
 * \file Color.cpp
 * 04/02/2010 KR Initial implementation
 *
 */


#ifndef COLOR_H_
#define COLOR_H_

#include "TemplatePoint.h"

#include <fstream>
#include <cstdlib>
#include <map>
#include <list>

namespace GEOLIB {

typedef TemplatePoint<unsigned char> Color;

/// Returns a random RGB colour.
Color* getRandomColor();

/// Reads a color-lookup-table from a file and writes it to a map.
int readColorLookupTable(std::map<std::string, GEOLIB::Color*> &colors, const std::string &filename);

/// Uses a color-lookup-table (in form of a map) to return a colour for a specified name. If the name is not
/// in the colortable a new entry is created with the new name and a random colour.
const Color* getColor(const std::string &id, std::map<std::string, GEOLIB::Color*> &colors);

/// Convenience function to use the getColor method with numbers as identifiers.
const Color* getColor(double id, std::map<std::string, GEOLIB::Color*> &colors);


} // namespace

#endif /* COLOR_H_ */
