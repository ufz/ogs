/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CoordinateSystem.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "BaseLib/Logging.h"

namespace FileIO
{
namespace Gocad
{
std::string parseName(std::string const& str)
{
    std::string name;
    std::size_t const start = str.find_first_of('\"');
    if (start != std::string::npos)
    {
        std::size_t const end = str.find_last_of('\"');
        name = str.substr(start + 1, end - start - 1);
    }
    else
    {
        name = str.substr(str.find_first_of(' '), str.length());
    }
    boost::algorithm::trim(name);
    return name;
}

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
    name = parseName(line);
    projection = "";
    datum = "";

    while (std::getline(in, line))
    {
        tok.assign(line);
        it = tok.begin();
        if (*it == "AXIS_NAME")
        {
            ++it;
            axis_name_u = *it;
            ++it;
            axis_name_v = *it;
            ++it;
            axis_name_w = *it;
        }
        else if (*it == "AXIS_UNIT")
        {
            ++it;
            axis_unit_u = *it;
            ++it;
            axis_unit_v = *it;
            ++it;
            axis_unit_w = *it;
        }
        else if (*it == "ZPOSITIVE")
        {
            ++it;
            if (*it == "Depth")
            {
                z_positive = ZPOSITIVE::Depth;
            }
            else
            {
                z_positive = ZPOSITIVE::Elevation;
            }
        }
        else if (*it == "PROJECTION")
        {
            projection = parseName(line);
        }
        else if (*it == "DATUM")
        {
            datum = parseName(line);
        }
        else if (*it == "END_ORIGINAL_COORDINATE_SYSTEM")
        {
            return true;
        }
        else
        {
            WARN("CoordinateSystem::parse() - Unknown keyword found: {:s}",
                 line);
        }
    }
    ERR("Error: Unexpected end of file.");
    return false;
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
