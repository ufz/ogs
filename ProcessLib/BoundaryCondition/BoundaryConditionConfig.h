/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "GeoLib/GEOObjects.h"

namespace ProcessLib
{

struct BoundaryConditionConfig final
{
    BoundaryConditionConfig(BaseLib::ConfigTree&& config_,
                            GeoLib::GeoObject const& geometry_,
                            int const component_id_)
        : config(std::move(config_)),
          geometry(geometry_),
          component_id(component_id_)
    {
    }

    BoundaryConditionConfig(BoundaryConditionConfig&& other)
        : config(std::move(other.config)),
          geometry(other.geometry),
          component_id(other.component_id)
    {
    }

    BaseLib::ConfigTree config;
    GeoLib::GeoObject const& geometry;
    int const component_id;
};

}  // ProcessLib
