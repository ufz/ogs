/*
 * GMSHPoint.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: TF
 */

#include <cmath>
#include <limits>

#include "MeshIO/GMSHPoint.h"

namespace FileIO {

GMSHPoint::GMSHPoint(GEOLIB::Point const& pnt, size_t id, double mesh_density) :
	GEOLIB::PointWithID(pnt, id), _mesh_density(mesh_density)
{}

void GMSHPoint::write(std::ostream &os) const
{
	os << "Point(" << _id << ") = {" << _x[0] << "," << _x[1] << ", 0.0";
	if (fabs(_mesh_density) > std::numeric_limits<double>::epsilon()) {
		os << ", " << _mesh_density << "};";
	} else {
		os << "};";
	}
}

GMSHPoint::~GMSHPoint()
{}

std::ostream& operator<< (std::ostream &os, GMSHPoint const& p)
{
	p.write (os);
	return os;
}

} // end namespace FileIO
