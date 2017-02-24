/**
 *  \brief A function for creating viscosity model
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file  CreateViscosityModel.cpp
 *
 */

#include "CreateViscosityModel.h"

#include "BaseLib/Error.h"

#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "LinearPressureDependentViscosity.h"
#include "TemperatureDependentViscosity.h"
#include "VogelsLiquidDynamicViscosity.h"
#include "WaterViscosityIAPWS.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 *     \param config  ConfigTree object which contains the input data
 *                    including `<type>LinearPressure</type>`
 *                    and it has a tag of `<viscosity>`
 */
static std::unique_ptr<FluidProperty> createLinearPressureDependentViscosity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__viscosity__type}
    config.checkConfigParameter("type", "LinearPressure");

    //! \ogs_file_param{material__fluid__viscosity__LinearPressure__mu0}
    const double mu0 = config.getConfigParameter<double>("mu0");

    //! \ogs_file_param{material__fluid__viscosity__LinearPressure__p0}
    const double p0 = config.getConfigParameter<double>("p0");

    //! \ogs_file_param{material__fluid__viscosity__LinearPressure__gamma}
    const double gamma = config.getConfigParameter<double>("gamma");

    return std::unique_ptr<FluidProperty>(
        new LinearPressureDependentViscosity(mu0, p0, gamma));
}

/**
 *     \param config  ConfigTree object which contains the input data
 *                    `<type>TemperatureDependent</type>`
 *                     and it has a tag of `<viscosity>`
 */
static std::unique_ptr<FluidProperty> createTemperatureDependentViscosity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__viscosity__type}
    config.checkConfigParameter("type", "TemperatureDependent");

    //! \ogs_file_param{material__fluid__viscosity__TemperatureDependent__mu0}
    const double mu0 = config.getConfigParameter<double>("mu0");

    //! \ogs_file_param{material__fluid__viscosity__TemperatureDependent__tc}
    const double Tc = config.getConfigParameter<double>("tc");

    //! \ogs_file_param{material__fluid__viscosity__TemperatureDependent__tv}
    const double Tv = config.getConfigParameter<double>("tv");

    return std::unique_ptr<FluidProperty>(
        new TemperatureDependentViscosity(mu0, Tc, Tv));
}

std::unique_ptr<FluidProperty> createViscosityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__viscosity__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        //! \ogs_file_param{material__fluid__viscosity__type}
        config.checkConfigParameter("type", "Constant");
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__viscosity__Constant__value}
            config.getConfigParameter<double>("value")));
    }
    else if (type == "LinearPressure")
        return createLinearPressureDependentViscosity(config);
    else if (type == "TemperatureDependent")
        return createTemperatureDependentViscosity(config);
    else if (type == "Vogels")
    {
        //! \ogs_file_param{material__fluid__viscosity__type}
        config.checkConfigParameter("type", "Vogels");

        INFO("Using Vogels model, which gives viscosity in SI unit, Pa s");
        auto const fluid_type =
            //! \ogs_file_param{material__fluid__viscosity__Vogels__liquid_type}
            config.peekConfigParameter<std::string>("liquid_type");
        if (fluid_type == "Water")
        {
            //! \ogs_file_param{material__fluid__viscosity__Vogels__liquid_type}
            config.checkConfigParameter("liquid_type", "Water");

            const VogelsViscosityConstantsWater constants;
            return std::unique_ptr<FluidProperty>(
                new VogelsLiquidDynamicViscosity<VogelsViscosityConstantsWater>(
                    constants));
        }
        else if (fluid_type == "CO2")
        {
            //! \ogs_file_param{material__fluid__viscosity__Vogels__liquid_type}
            config.checkConfigParameter("liquid_type", "CO2");
            const VogelsViscosityConstantsCO2 constants;
            return std::unique_ptr<FluidProperty>(
                new VogelsLiquidDynamicViscosity<VogelsViscosityConstantsCO2>(
                    constants));
        }
        else if (fluid_type == "CH4")
        {
            //! \ogs_file_param{material__fluid__viscosity__Vogels__liquid_type}
            config.checkConfigParameter("liquid_type", "CH4");
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
    else if (type == "WaterViscosityIAPWS")
    {
        //! \ogs_file_param{material__fluid__viscosity__type}
        config.checkConfigParameter("type", "WaterViscosityIAPWS");
        return std::unique_ptr<FluidProperty>(new WaterViscosityIAPWS());
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
