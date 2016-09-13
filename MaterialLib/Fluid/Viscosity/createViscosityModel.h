/*!
   \file  createViscosityModel.h

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
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
/// Create a viscosity model
/// \param config  ConfigTree object has a tag of <viscosity>
std::unique_ptr<FluidProperty> createViscosityModel(BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace
#endif
