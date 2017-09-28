/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "cWater.h"
#include "../Properties/properties.h"

namespace MaterialPropertyLib
{
Water::Water()
{
    _properties[name] = std::make_unique<Constant>("Water");
};

}  // MaterialPropertyLib
