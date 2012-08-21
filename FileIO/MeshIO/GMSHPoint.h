/*
 * GMSHPoint.h
 *
 *  Created on: Mar 21, 2012
 *      Author: TF
 */

#ifndef GMSHPOINT_H_
#define GMSHPOINT_H_

// GEOLIB
#include "PointWithID.h"

namespace FileIO {

class GMSHPoint : public GeoLib::PointWithID {
public:
	GMSHPoint(GeoLib::Point const& pnt, size_t id, double mesh_density);
	virtual ~GMSHPoint();
	void write(std::ostream &os) const;
private:
	double _mesh_density;
};

/** overload the output operator for class GMSHPoint */
std::ostream& operator<< (std::ostream &os, GMSHPoint const& p);

}

#endif /* GMSHPOINT_H_ */
