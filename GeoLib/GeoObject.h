/*
 * GeoObject.h
 *
 *  Created on: Aug 27, 2010
 *      Author: TF
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GEOLIB {

/**
 * \ingroup GEOLIB
 *
 * \brief Base class for classes Point, Polyline, Surface.
 */

class GeoObject {
public:
	GeoObject() {};
	virtual ~GeoObject() {};
};

} // end namespace GEOLIB

#endif /* GEOOBJECT_H_ */
