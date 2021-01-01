/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cstddef>
#include <string>

#include "BaseLib/Logging.h"

namespace FileIO
{
namespace Gocad
{
struct Property final
{
    std::size_t _property_id{};
    std::string _property_name;
    std::string _property_class_name;
    std::string _property_unit;
    std::string _property_data_type;
    std::string _property_data_fname;
    double _property_no_data_value{};

    bool checkID(std::string const& id_string) const
    {
        if (_property_id != std::stoul(id_string))
        {
            ERR("Expected property id {:d} but found {:d}.",
                _property_id,
                std::stoul(id_string));
            return false;
        }
        return true;
    }

    std::vector<double> _property_data;
};

std::ostream& operator<<(std::ostream& os, Property const& p);

Property parseGocadPropertyMetaData(std::string& line, std::istream& in,
                                    std::string const& path);
}  // end namespace Gocad
}  // end namespace FileIO
