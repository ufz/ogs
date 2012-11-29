/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GMSHPoint.h
 *
 *  Created on 2012-03-21 by Thomas Fischer
 */

#ifndef GMSHPOINT_H_
#define GMSHPOINT_H_

// GeoLib
#include "PointWithID.h"

namespace FileIO {

class GMSHPoint : public GeoLib::PointWithID {
public:
	GMSHPoint(GeoLib::Point const& pnt, std::size_t id, double mesh_density);
	virtual ~GMSHPoint();
	void write(std::ostream &os) const;
private:
	double _mesh_density;
};

/** overload the output operator for class GMSHPoint */
std::ostream& operator<< (std::ostream &os, GMSHPoint const& p);

}

#endif /* GMSHPOINT_H_ */
