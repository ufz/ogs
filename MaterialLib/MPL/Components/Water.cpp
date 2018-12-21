/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Water.h"
#include "MaterialLib/MPL/Properties/Properties.h"

namespace MaterialPropertyLib
{
Water::Water(std::unique_ptr<PropertyArray>&& properties)
{
    _properties[PropertyType::name] = std::make_unique<Constant>("Water");

    if (properties)
    {
        overwriteExistingProperties(_properties, *properties);
    }
}
}  // namespace MaterialPropertyLib
