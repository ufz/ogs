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
/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>fluid</type> and it has
                   a tag of <viscosity>
*/
static std::unique_ptr<FluidProperty> createLinearPressureDependentViscosity(
    BaseLib::ConfigTree const& config)
{
    std::array<double, 3> parameters = {
        {//! \ogs_file_param{material__fluid__viscosity__temperature_dependent__mu0}
         config.getConfigParameter<double>("mu0"),
         //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__p0}
         config.getConfigParameter<double>("p0"),
         //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__gamma}
         config.getConfigParameter<double>("gamma")}};
    return std::unique_ptr<FluidProperty>(
        new LinearPressureDependentViscosity(parameters));
}

/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>fluid</type> and it has
                   a tag of <viscosity>
*/
static std::unique_ptr<FluidProperty> createTemperatureDependentViscosity(
    BaseLib::ConfigTree const& config)
{
    std::array<double, 3> parameters = {
        {//! \ogs_file_param{material__fluid__viscosity__temperature_dependent__mu0}
         config.getConfigParameter<double>("mu0"),
         //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tc}
         config.getConfigParameter<double>("tc"),
         //! \ogs_file_param{material__fluid__viscosity__temperature_dependent__tv}
         config.getConfigParameter<double>("tv")}};
    return std::unique_ptr<FluidProperty>(
        new TemperatureDependentViscosity(parameters));
}

std::unique_ptr<FluidProperty> createViscosityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__viscosity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__viscosity__Constant__value}
            config.getConfigParameter<double>("value")));
    //! \ogs_file_param{material__fluid__viscosity__LinearPressure}
    else if (type == "LinearPressure")
        return createLinearPressureDependentViscosity(config);
    //! \ogs_file_param{material__fluid__viscosity__TemperatureDependent}
    else if (type == "TemperatureDependent")
        return createTemperatureDependentViscosity(config);
    //! \ogs_file_param{material__fluid__viscosity__Vogels}
    else if (type == "Vogels")
    {
        auto const fluid_type =
            //! \ogs_file_param{material__fluid__viscosity__Vogels__fluid_type}
            config.getConfigParameter<std::string>("liquid_type");
        //! \ogs_file_param{material__fluid__viscosity__Vogels__Water}
        if (fluid_type == "Water")
        {
            const VogelsViscosityConstantsWater constants;
            return std::unique_ptr<FluidProperty>(
                new VogelsLiquidDynamicViscosity<VogelsViscosityConstantsWater>(
                    constants));
        }
        //! \ogs_file_param{material__fluid__viscosity__Vogels__CO2}
        else if (fluid_type == "CO2")
        {
            const VogelsViscosityConstantsCO2 constants;
            return std::unique_ptr<FluidProperty>(
                new VogelsLiquidDynamicViscosity<VogelsViscosityConstantsCO2>(
                    constants));
        }
        //! \ogs_file_param{material__fluid__viscosity__Vogels__CH4}
        else if (fluid_type == "CH4")
        {
            const VogelsViscosityConstantsCH4 constants;
            return std::unique_ptr<FluidProperty>(
                new VogelsLiquidDynamicViscosity<VogelsViscosityConstantsCH4>(
                    constants));
        }
        else
        {
            OGS_FATAL(
                "The fluid type %s for Vogels model is unavailable.\n"
                "The available fluid types are Water, CO2 and CH4\n",
                fluid_type.data());
        }
    }
    else
    {
        OGS_FATAL(
            "The viscosity type %s is unavailable.\n"
            "The available types are \n\tConstant, \n\tLinearPressure "
            "\n\tTemperatureDependent, \n\tVogels\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
