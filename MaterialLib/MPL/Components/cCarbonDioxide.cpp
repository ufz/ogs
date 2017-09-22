/**
 * \author Norbert Grunwald
 * \date   12.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "cCarbonDioxide.h"
#include "../Properties/properties.h"

namespace MaterialPropertyLib
{

CarbonDioxide::CarbonDioxide()
{
	_properties[name] = std::make_unique<Constant>("CarbonDioxide");
};

} // MaterialPropertyLib
