/**
 * \file
 * \author Thomas Fischer
 * \date   2010-08-27
 * \brief  Base class for classes Point, Polyline, Surface.
 * \ingroup GeoLib
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GeoLib
{
class GeoObject
{
public:
	GeoObject() {}
	virtual ~GeoObject() {}
};
} // end namespace GeoLib

#endif /* GEOOBJECT_H_ */
