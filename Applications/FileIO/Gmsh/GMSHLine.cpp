// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
