/**
 * \file
 * \author Thomas Fischer
 * \date   2010-08-27
 * \brief  Base class for classes Point, Polyline, Surface.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "GeoType.h"

namespace GeoLib
{
struct GeoObject
{
    virtual ~GeoObject() = default;
    /// return a geometry type
    virtual GEOTYPE getGeoType() const = 0;
};
}  // end namespace GeoLib
