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
#include "GMSHLineLoop.h"

namespace FileIO
{
namespace GMSH
{
GMSHLineLoop::GMSHLineLoop(bool is_sfc) : is_sfc_(is_sfc) {}

GMSHLineLoop::~GMSHLineLoop()
{
    const std::size_t n_lines(lines_.size());
    for (std::size_t k(0); k < n_lines; k++)
    {
        delete lines_[k];
    }
}

void GMSHLineLoop::write(std::ostream& os, std::size_t line_offset,
                         std::size_t sfc_offset) const
{
    const std::size_t n_lines(lines_.size());
    for (std::size_t k(0); k < n_lines; k++)
    {
        (lines_[k])->write(os, line_offset + k);
    }
    os << "Line Loop(" << line_offset + n_lines << ") = {";
    for (std::size_t k(0); k < n_lines - 1; k++)
    {
        os << line_offset + k << ",";
    }
    os << line_offset + n_lines - 1 << "};\n";

    if (is_sfc_)
    {
        // write plane surface
        os << "Plane Surface (" << sfc_offset << ") = {"
           << line_offset + n_lines << "};\n";
    }
}

}  // end namespace GMSH
}  // end namespace FileIO
