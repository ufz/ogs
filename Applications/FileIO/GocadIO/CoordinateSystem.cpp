/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <iostream>

#include <boost/tokenizer.hpp>

#include "CoordinateSystem.h"

namespace FileIO
{
namespace Gocad
{
bool CoordinateSystem::parse(std::istream& in)
{
    std::string line;  // string used for reading a line

    boost::char_separator<char> sep("-;| \"");

    std::getline(in, line);  // NAME name
    boost::tokenizer<boost::char_separator<char>> tok(line, sep);
    auto it(tok.begin());
    if (*it != "NAME")
    {
        return false;
    }
    it++;
    name = *it;

    std::getline(in, line);  // AXIS_NAME "n1" "n2" "n3"
    tok.assign(line);
    it = tok.begin();
    if (*it != "AXIS_NAME")
    {
        return false;
    }
    ++it;
    axis_name_u = *it;
    ++it;
    axis_name_v = *it;
    ++it;
    axis_name_w = *it;

    std::getline(in, line);  // AXIS_UNIT "u1" "u2" "u3"
    tok.assign(line);
    it = tok.begin();
    if (*it != "AXIS_UNIT")
    {
        return false;
    }
    ++it;
    axis_unit_u = *it;
    ++it;
    axis_unit_v = *it;
    ++it;
    axis_unit_w = *it;

    std::getline(in, line);  // ZPOSITIVE Depth || Elevation
    tok.assign(line);
    it = tok.begin();
    if (*it != "ZPOSITIVE")
    {
        return false;
    }
    ++it;
    if (*it != "Depth")
    {
        z_positive = ZPOSITIVE::Depth;
    }
    else
    {
        z_positive = ZPOSITIVE::Elevation;
    }

    std::getline(in, line);  // END_ORIGINAL_COORDINATE_SYSTEM
    tok.assign(line);
    it = tok.begin();
    return *it == "END_ORIGINAL_COORDINATE_SYSTEM";
}

std::ostream& operator<<(std::ostream& os, CoordinateSystem const& c)
{
    os << "Gocad CoordinateSystem " << c.name
       << "\nAxis names: " << c.axis_name_u << ", " << c.axis_name_v << ", "
       << c.axis_name_w << "\nAxis units: " << c.axis_unit_u << ", "
       << c.axis_unit_v << ", " << c.axis_unit_w << "\nZ-orientation: ";
    if (c.z_positive == FileIO::Gocad::CoordinateSystem::ZPOSITIVE::Depth)
    {
        os << "downwards";
    }
    else
    {
        os << "upwards";
    }
    os << "\n";
    return os;
}

}  // end namespace Gocad
}  // end namespace FileIO
