/*!
   \file  createFluidDensityModel.h
   \brief create an instance of a fluid density class.

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef CREATE_FLUID_DENSITY_MODEL_H_
#define CREATE_FLUID_DENSITY_MODEL_H_

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/// Create a density model
/// \param config  ConfigTree object has a tag of <density>
FluidProperty* createFluidDensityModel(BaseLib::ConfigTree const* const config);
}
}  // end namespace
#endif
