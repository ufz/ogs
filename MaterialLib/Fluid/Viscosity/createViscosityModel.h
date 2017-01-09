/**
 *  \brief A function for creating viscosity model
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file  createViscosityModel.h
 *
 */

#ifndef CREATE_VISCOSITY_MODEL_H_
#define CREATE_VISCOSITY_MODEL_H_

#include <memory>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
///  Create a viscosity model
///  \param config  ConfigTree object has a tag of `<viscosity>`
std::unique_ptr<FluidProperty> createViscosityModel(
    BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace
#endif
