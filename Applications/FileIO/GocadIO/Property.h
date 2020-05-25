/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::size_t property_id_{};
    std::string property_name_;
    std::string property_class_name_;
    std::string property_unit_;
    std::string property_data_type_;
    std::string property_data_fname_;
    double property_no_data_value_{};

    bool checkID(std::string const& id_string) const
    {
        if (property_id_ != std::stoul(id_string))
        {
            ERR("Expected property id {:d} but found {:d}.",
                property_id_,
                std::stoul(id_string));
            return false;
        }
        return true;
    }

    std::vector<double> property_data_;
};

std::ostream& operator<<(std::ostream& os, Property const& p);

Property parseGocadPropertyMetaData(std::string& line, std::istream& in,
                                    std::string const& path);
}  // end namespace Gocad
}  // end namespace FileIO
