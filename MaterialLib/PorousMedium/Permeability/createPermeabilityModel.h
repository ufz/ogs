/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createPermeabilityModel.h
 *
 * Created on August 17, 2016, 2:36 PM
 */

#pragma once

#include <memory>

#include "Permeability.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
/** Create a porosity model
 * @param config  ConfigTree object has a tag of `<permeability>` that
 * describes the permeability relationsship and contains the name of the
 * parameter
 * @param parameters a vector containing the available parameters
 */
std::unique_ptr<Permeability> createPermeabilityModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);

}  // end of namespace
}  // end of namespace
