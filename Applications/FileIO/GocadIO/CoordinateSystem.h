// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iosfwd>
#include <string>

namespace FileIO
{
namespace Gocad
{

class CoordinateSystem final
{
public:
    bool parse(std::istream & in);

    enum class ZPOSITIVE
    {
        Depth,     // z is increasing downwards
        Elevation  // z is increasing upwards
    };

    std::string name;
    std::string projection;
    std::string datum;
    std::string axis_name_u, axis_name_v, axis_name_w;
    std::string axis_unit_u, axis_unit_v, axis_unit_w;
    ZPOSITIVE z_positive = ZPOSITIVE::Elevation;
};

std::ostream& operator<<(std::ostream& os, CoordinateSystem const& c);

}  // end namespace Gocad
}  // end namespace FileIO
