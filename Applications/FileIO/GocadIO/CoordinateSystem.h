/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <fstream>

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
    std::string axis_name_u, axis_name_v, axis_name_w;
    std::string axis_unit_u, axis_unit_v, axis_unit_w;
    ZPOSITIVE z_positive;
};

std::ostream& operator<<(std::ostream& os, CoordinateSystem const& c);

}  // end namespace Gocad
}  // end namespace FileIO
