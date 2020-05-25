/**
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>
#include "GMSHLine.h"

namespace FileIO
{
namespace GMSH
{

GMSHLine::GMSHLine(std::size_t start_point_id, std::size_t end_point_id) :
    start_pnt_id_(start_point_id), end_pnt_id_(end_point_id)
{}

void GMSHLine::write(std::ostream &os, std::size_t id) const
{
    os << "Line(" << id << ") = {" << start_pnt_id_ << "," << end_pnt_id_ << "};\n";
}

} // end namespace GMSH
} // end namespace FileIO
