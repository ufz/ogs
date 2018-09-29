/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

GMSHLineLoop::GMSHLineLoop(bool is_sfc) :
    _is_sfc(is_sfc)
{}

GMSHLineLoop::~GMSHLineLoop()
{
    const std::size_t n_lines (_lines.size());
    for (std::size_t k(0); k<n_lines; k++) {
        delete _lines[k];
    }
}

void GMSHLineLoop::addLine(GMSHLine* line)
{
    _lines.push_back(line);
}

void GMSHLineLoop::write(std::ostream &os, std::size_t line_offset, std::size_t sfc_offset) const
{
    const std::size_t n_lines (_lines.size());
    for (std::size_t k(0); k<n_lines; k++) {
        (_lines[k])->write(os, line_offset+k);
    }
    os << "Line Loop(" << line_offset+n_lines << ") = {";
    for (std::size_t k(0); k < n_lines - 1; k++)
        os << line_offset + k << ",";
    os << line_offset + n_lines - 1 << "};\n";

    if (_is_sfc) {
        // write plane surface
        os << "Plane Surface (" << sfc_offset << ") = {" << line_offset+n_lines << "};\n";
    }

}

}  // end namespace GMSH
}  // end namespace FileIO
