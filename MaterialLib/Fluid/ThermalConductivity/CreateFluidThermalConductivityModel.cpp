/**
 * \file
 *  \brief A function for creating a thermal conductivity model for fluid
 *
 *  \copyright
 *   Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateFluidThermalConductivityModel.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__thermal_conductivity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        return std::make_unique<ConstantFluidProperty>(
            //! \ogs_file_param{material__fluid__thermal_conductivity__Constant__value}
            config.getConfigParameter<double>("value"));
    }
    // TODO: add more models

    OGS_FATAL(
        "The viscosity type {:s} is unavailable.\n"
        "The available type is \n\tConstant\n",
        type.data());
}

}  // namespace Fluid
}  // namespace MaterialLib
