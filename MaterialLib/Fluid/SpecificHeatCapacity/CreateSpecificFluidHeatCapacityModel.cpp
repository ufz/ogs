/**
 *  \brief A function for creating a specific heat capacity model for fluid
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file CreateSpecificFluidHeatCapacityModel.cpp
 *
 */

#include "CreateSpecificFluidHeatCapacityModel.h"

#include "BaseLib/Error.h"
#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
std::unique_ptr<FluidProperty> createSpecificFluidHeatCapacityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__specific_heat_capacity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::make_unique<ConstantFluidProperty>(
            //! \ogs_file_param{material__fluid__specific_heat_capacity__Constant__value}
            config.getConfigParameter<double>("value"));
    // TODO: add more models

    OGS_FATAL(
        "The viscosity type %s is unavailable.\n"
        "The available type is \n\tConstant\n",
        type.data());
}

}  // end namespace
}  // end namespace
