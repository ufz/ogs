/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 16, 2016, 1:16 PM
 */

#pragma once

#include <memory>
#include "ParameterLib/Parameter.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;

/** Create a porosity model
 *  @param config  ConfigTree object has a tag of `<porosity>` that describes
 *  the porosity relationsship and contains the name of the parameter
 *  @param parameters a vector containing the available parameters
 */
std::unique_ptr<Porosity> createPorosityModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace PorousMedium
}  // namespace MaterialLib
