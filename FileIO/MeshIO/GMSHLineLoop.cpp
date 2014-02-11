/**
 * \file
 * \author Thomas Fischer
 * \date   Mar 22, 2012
 * \brief  Implementation of the GMSHLineLoop class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshIO/GMSHLineLoop.h"

namespace FileIO 
{
namespace GMSH {

GMSHLineLoop::GMSHLineLoop(bool is_sfc) :
	_is_sfc(is_sfc)
{}

GMSHLineLoop::~GMSHLineLoop()
{
	const size_t n_lines (_lines.size());
	for (size_t k(0); k<n_lines; k++) {
		delete _lines[k];
	}
}

void GMSHLineLoop::addLine(GMSHLine* line)
{
	_lines.push_back(line);
}

void GMSHLineLoop::write(std::ostream &os, size_t line_offset, size_t sfc_offset) const
{
	const size_t n_lines (_lines.size());
	for (size_t k(0); k<n_lines; k++) {
		(_lines[k])->write(os, line_offset+k);
	}
	os << "Line Loop(" << line_offset+n_lines << ") = {";
	for (size_t k(0); k < n_lines - 1; k++)
		os << line_offset + k << ",";
	os << line_offset + n_lines - 1 << "};\n";

	if (_is_sfc) {
		// write plane surface
		os << "Plane Surface (" << sfc_offset << ") = {" << line_offset+n_lines << "};\n";
	}

}

}
} // end namespace FileIO
