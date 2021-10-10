/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Property.h"

#include <algorithm>
#include <boost/tokenizer.hpp>
#include <iterator>
#include <sstream>

#include "BaseLib/StringTools.h"

namespace FileIO
{
namespace Gocad
{
std::ostream& operator<<(std::ostream& os, Property const& p)
{
    return os << "'" << p._property_name << "' '" << p._property_id << "' '"
              << p._property_data_type << "' '" << p._property_data_fname
              << "'";
}

Property parseGocadPropertyMetaData(std::string& line, std::istream& in,
                                    std::string const& path)
{
    boost::char_separator<char> sep("\t ");
    boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
    auto tok_it(tokens.begin());
    // A property section in Gocad file starts with a line
    // PROPERTY id "property_name"
    if (*tok_it != "PROPERTY")
    {
        ERR("Expected PROPERTY keyword but '{:s}' found.", tok_it->c_str());
        throw std::runtime_error(
            "In parseGocadPropertyMetaData() expected PROPERTY keyword not "
            "found.");
    }
    tok_it++;

    Property prop;
    prop._property_id = std::stoul(*tok_it);
    tok_it++;
    prop._property_name = *tok_it;
    tok_it++;
    while (tok_it != tokens.end())
    {
        prop._property_name += " " + *tok_it;
        tok_it++;
    }
    BaseLib::trim(prop._property_name, '\"');

    auto checkPropertyID =
        [](boost::tokenizer<boost::char_separator<char>>::iterator const&
               tok_it,
           Property const& prop)
    {
        if (!prop.checkID(*tok_it))
        {
            throw std::runtime_error(
                "parseGocadPropertyMetaData(): id mismatch.");
        }
    };

    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (line.back() == '\r')
        {
            line.pop_back();
        }
        tokens.assign(line);

        tok_it = tokens.begin();
        // this is the last entry of the property
        if (*tok_it == "PROP_FILE")
        {
            checkPropertyID(++tok_it, prop);
            tok_it++;
            std::string tmp(*tok_it);
            tok_it++;
            while (tok_it != tokens.end())
            {
                tmp += " " + *tok_it;
                tok_it++;
            }
            BaseLib::trim(tmp, '\"');
            if (tmp.front() == ' ')
            {
                tmp.erase(tmp.begin());
            }
            prop._property_data_fname = path + tmp;
            return prop;
        }

        if (*tok_it == "PROPERTY_CLASS")
        {
            checkPropertyID(++tok_it, prop);
            tok_it++;
            prop._property_class_name = *tok_it;
        }

        if (*tok_it == "PROPERTY_SUBCLASS")
        {
            checkPropertyID(++tok_it, prop);
            tok_it++;
            if (*tok_it != "QUANTITY" && *tok_it != "ENUM")
            {
                ERR("Expected keywords QUANTITY or ENUM, but found '{:s}'.",
                    tok_it->c_str());
                throw std::runtime_error(
                    "parseGocadPropertyMetaData(): Expected keywords QUANTITY "
                    "or ENUM, but found '" +
                    *tok_it + "'.");
            }
            if (*tok_it == "QUANTITY")
            {
                tok_it++;
                prop._property_data_type = *tok_it;
            }
        }

        if (*tok_it == "PROP_UNIT" || *tok_it == "PROP_ORIGINAL_UNIT")
        {
            checkPropertyID(++tok_it, prop);
            tok_it++;
            prop._property_unit = *tok_it;
        }

        if (*tok_it == "PROP_NO_DATA_VALUE")
        {
            checkPropertyID(++tok_it, prop);
            tok_it++;
            prop._property_no_data_value = std::stoul(*tok_it);
        }
    }
    return prop;
}

}  // end namespace Gocad
}  // end namespace FileIO
