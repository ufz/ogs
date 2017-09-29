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
#pragma once

#include "../mpComponent.h"

namespace MaterialPropertyLib
{
/**
 * \class Water
 * \brief A class for Water derived from Component
 * \details This class can holds material constants and default
 * properties of ordinary water
 */
class Water final : public Component
{
public:
    Water();
};

}  // MaterialPropertyLib
