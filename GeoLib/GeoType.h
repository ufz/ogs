/*
 * GeoType.h
 *
 *  Created on: Jun 17, 2010
 *      Author: TF
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_

#include <string>

namespace GeoLib {

/**
 * \ingroup GeoLib
 */

enum GEOTYPE {
	INVALID = 0,
	POINT,     //!< POINT
	POLYLINE,  //!< POLYLINE
	SURFACE,   //!< SURFACE
	VOLUME,    //!< VOLUME
	GEODOMAIN, //!< GEODOMAIN
	COLUMN     //!< COLUMN. //WW/JOD 
};

GEOTYPE convertGeoType (const std::string& geo_type_str);

std::string convertGeoTypeToString (GEOTYPE geo_type);

} // end namespace GeoLib

#endif /* GEOTYPE_H_ */
