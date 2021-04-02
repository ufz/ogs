/**
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GMSHLine.h"

#include <ostream>

namespace FileIO
{
namespace GMSH
{
GMSHLine::GMSHLine(std::size_t start_point_id, std::size_t end_point_id)
    : _start_pnt_id(start_point_id), _end_pnt_id(end_point_id)
{
}

void GMSHLine::write(std::ostream& os, std::size_t id) const
{
    os << "Line(" << id << ") = {" << _start_pnt_id << "," << _end_pnt_id
       << "};\n";
}

}  // end namespace GMSH
}  // end namespace FileIO
