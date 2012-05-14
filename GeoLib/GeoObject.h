/*
 * GeoObject.h
 *
 *  Created on: Aug 27, 2010
 *      Author: TF
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GeoLib {

/**
 * \ingroup GeoLib
 *
 * \brief Base class for classes Point, Polyline, Surface.
 */

class GeoObject {
public:
	GeoObject() {};
	virtual ~GeoObject() {};
};

} // end namespace GeoLib

#endif /* GEOOBJECT_H_ */
