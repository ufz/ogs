/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TINInterface.h"

#include <fstream>
#include <limits>

#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/Surface.h"


namespace FileIO
{

GeoLib::Surface* TINInterface::readTIN(std::string const& fname, std::vector<GeoLib::Point*> &pnt_vec, std::vector<std::string>* errors)
{
	// open file
	std::ifstream in(fname.c_str());
	if (!in) {
		WARN("readTIN(): could not open stream from %s.", fname.c_str());
		if (errors) errors->push_back ("readTINFile error opening stream from " + fname);
		return nullptr;
	}

	GeoLib::Surface* sfc = new GeoLib::Surface(pnt_vec);
	std::size_t id;
	double x, y, z;
	while (in)
	{
		// read id
		in >> id;
		// determine size
		std::size_t pnt_pos(pnt_vec.size());
		// read first point
		in >> x >> y >> z;
		pnt_vec.push_back(new GeoLib::Point(x, y, z));
		// read second point
		in >> x >> y >> z;
		pnt_vec.push_back(new GeoLib::Point(x, y, z));
		// read third point
		in >> x >> y >> z;
		pnt_vec.push_back(new GeoLib::Point(x, y, z));
		// create new Triangle
		sfc->addTriangle(pnt_pos, pnt_pos + 1, pnt_pos + 2);
	}

	if (sfc->getNTriangles() == 0) {
		WARN("readTIN(): No triangle found.", fname.c_str());
		if (errors) errors->push_back ("readTIN error because of no triangle found");
		delete sfc;
		return nullptr;
	}

	return sfc;
}

void TINInterface::writeSurfaceAsTIN(GeoLib::Surface const& surface, std::string const& file_name)
{
	std::ofstream os (file_name.c_str());
	os.precision(std::numeric_limits<double>::digits10);
	const std::size_t n_tris (surface.getNTriangles());
	for (std::size_t l(0); l < n_tris; l++) {
		GeoLib::Triangle const& tri (*(surface[l]));
		os << l << " " << *(tri.getPoint(0)) << " " << *(tri.getPoint(1)) << " " << *(tri.getPoint(2)) << "\n";
	}
	os.close();
}

} // end namespace GeoLib
