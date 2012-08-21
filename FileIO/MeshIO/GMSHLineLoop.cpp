/*
 * GMSHLineLoop.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: fischeth
 */

#include "MeshIO/GMSHLineLoop.h"

namespace FileIO {

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
	os << line_offset + n_lines - 1 << "};" << std::endl;

	if (_is_sfc) {
		// write plane surface
		os << "Plane Surface (" << sfc_offset << ") = {" << line_offset+n_lines << "};" << std::endl;
	}

}

} // end namespace FileIO
