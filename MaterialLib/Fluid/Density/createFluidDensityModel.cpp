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
FluidProperty* createFluidDensityModel(BaseLib::ConfigTree const* const config)
{
    //! \ogs_file_param{material__fluid__density__type}
    auto const type = config->getConfigParameter<std::string>("type");

    if (type.find("constant") != std::string::npos)
    {
        //! \ogs_file_param{material__fluid__density__value}
        return new ConstantFluidProperty(
            config->getConfigParameter<double>("value"));
    }
    else if (type.find("liquid_density") != std::string::npos)
        return new LiquidDensity(config);
    else if (type.find("temperature_dependent") != std::string::npos)
        return new LinearTemperatureDependentDensity(config);
    else if (type.find("ideal_gas_law") != std::string::npos)
    {
        //! \ogs_file_param{material__fluid__density__ideal_gas_law__malar_mass}
        return new IdealGasLaw(
            config->getConfigParameter<double>("molar_mass"));
    }
    else
    {
        OGS_FATAL(
            "The density type is unavailable.\n"
            "The available types are liquid_density "
            "temperature_dependent and ideal_gas_law");
    }
}

}  // end namespace
}  // end namespace
