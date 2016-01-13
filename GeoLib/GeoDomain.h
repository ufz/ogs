/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEODOMAIN_H_
#define GEODOMAIN_H_

#include "GeoType.h"
#include "GeoObject.h"

namespace GeoLib
{
class GeoDomain : public GeoObject
{
public:
	virtual ~GeoDomain() = default;
	/// return a geometry type
	virtual GEOTYPE getGeoType() const {return GEOTYPE::GEODOMAIN; }
};
} // end namespace GeoLib

#endif /* GEODOMAIN_H_ */
