/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "cWater.h"
#include "MaterialLib/MPL/Properties/properties.h"

namespace MaterialPropertyLib
{
Water::Water()
{
    _properties[name] = std::make_unique<Constant>("Water");
}
}  // namespace MaterialPropertyLib
