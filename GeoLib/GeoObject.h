/**
 * \file
 * \author Thomas Fischer
 * \date   2010-08-27
 * \brief  Base class for classes Point, Polyline, Surface.
 * \ingroup GeoLib
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

#include "GeoType.h"

namespace GeoLib
{
class GeoObject
{
public:
    GeoObject() = default;
    GeoObject(GeoObject const&) = default;
    GeoObject& operator=(GeoObject const&) = default;
    virtual ~GeoObject() = default;
    /// return a geometry type
    virtual GEOTYPE getGeoType() const = 0;
};
} // end namespace GeoLib

#endif /* GEOOBJECT_H_ */
