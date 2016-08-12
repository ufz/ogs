/*!
   \file  createViscosityModel.cpp

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "createViscosityModel.h"

#include "BaseLib/Error.h"

#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "LinearPressureDependentViscosity.h"
#include "TemperatureDependentViscosity.h"
#include "VogelsLiquidDynamicViscosity.h"

namespace MaterialLib
{
namespace Fluid
{
FluidProperty* createViscosityModel(BaseLib::ConfigTree const* const config)
{
    //! \ogs_file_param{material__fluid__viscosity__type}
    auto const type = config->getConfigParameter<std::string>("type");

    if (type.find("constant") != std::string::npos)
        //! \ogs_file_param{material__fluid__viscosity__value}
        return new ConstantFluidProperty(
            config->getConfigParameter<double>("value"));
    else if (type.find("linear_pressure") != std::string::npos)
        return new LinearPressureDependentViscosity(config);
    else if (type.find("temperature_dependent") != std::string::npos)
        return new TemperatureDependentViscosity(config);
    else if (type.find("vogels") != std::string::npos)
    {
        //! \ogs_file_param{material__fluid__viscosity__vogels__fluid_type}
        auto const fluid_type =
            config->getConfigParameter<std::string>("liquid_type");
        int type_id = -1;
        if (fluid_type.find("water") != std::string::npos)
            type_id = 0;
        else if (fluid_type.find("CO2") != std::string::npos)
            type_id = 1;
        else if (fluid_type.find("CH4") != std::string::npos)
            type_id = 2;

        if (type_id > -1 && type_id < 3)
        {
            return new VogelsLiquidDynamicViscosity(type_id);
        }
        else
        {
            OGS_FATAL(
                "The fluid type for Vogels model is unavailable.\n"
                "The available fluid types are water, CO2 and CH4 ");
        }
    }
    else
    {
        OGS_FATAL(
            "The viscosity type is unavailable.\n"
            "The available types are linear_pressure "
            "temperature_dependent and vogels");
    }
}

}  // end namespace
}  // end namespace
