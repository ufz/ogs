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
#ifndef MATERIALLIB_MPL_MPPROPERTY_H_
#define MATERIALLIB_MPL_MPPROPERTY_H_

#include "BaseLib/ConfigTree.h"
#include "mpEnums.h"
#include <array>

namespace MaterialPropertyLib
{
class Property
{
public:
    Property();
};

using PropertyArray = std::array<Property*, number_of_property_enums>;
Property* newProperty(BaseLib::ConfigTree const&);

} //MaterialPropertyLib


#endif /* MATERIALLIB_MPL_MPPROPERTY_H_ */
