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
#pragma once

#include "../mpComponent.h"

namespace MaterialPropertyLib
{
/**
 * \class CarbonDioxide
 * \brief A class for CO2 derived from Component
 * \details This class can holds material constants and default
 * properties of ordinary carbon dioxide.
 */
class CarbonDioxide : public Component
{
public:
    CarbonDioxide();
};

}  // MaterialPropertyLib
