/*!
   \file  createFluidDensityModel.cpp
   \brief create an instance of a fluid density class.

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "createFluidDensityModel.h"

#include "BaseLib/Error.h"

#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "IdealGasLaw.h"
#include "LinearTemperatureDependentDensity.h"
#include "LiquidDensity.h"

namespace MaterialLib
{
namespace Fluid
{
std::unique_ptr<FluidProperty> createFluidDensityModel(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__density__Constant_value}
            config.getConfigParameter<double>("value")) );
    }
    //! \ogs_file_param{material__fluid__density__LiquidDensity}
    else if (type == "LiquidDensity")
        return std::unique_ptr<FluidProperty>(new LiquidDensity(config));
    //! \ogs_file_param{material__fluid__density__TemperatureDependent}
    else if (type == "TemperatureDependent")
        return std::unique_ptr<FluidProperty>(
                new LinearTemperatureDependentDensity(config));
    //! \ogs_file_param{material__fluid__density__IdealGasLaw}
    else if (type == "IdealGasLaw")
    {
        //! \ogs_file_param{material__fluid__density__IdealGasLaw__molar_mass}
        return std::unique_ptr<FluidProperty>(new IdealGasLaw(
            config.getConfigParameter<double>("molar_mass")) );
    }
    else
    {
        OGS_FATAL(
            "The density type %s is unavailable.\n", type.data(),
            "The available types are Constant, LiquidDensity, "
            "TemperatureDependent and IdealGasLaw");
    }
}

}  // end namespace
}  // end namespace
