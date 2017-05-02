/*!
   \file  CreateFluidDensityModel.cpp
   \brief create an instance of a fluid density class.

   \copyright
    Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <array>

#include "CreateFluidDensityModel.h"

#include "BaseLib/Error.h"

#include "IdealGasLaw.h"
#include "LinearConcentrationDependentDensity.h"
#include "LinearTemperatureDependentDensity.h"
#include "LiquidDensity.h"
#include "WaterDensityIAPWSIF97Region1.h"

#include "MaterialLib/Fluid/ConstantFluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 *  \param config  ConfigTree object which contains the input data
 *                 including `<type>LiquidDensity</type>`
 *                 and it has a tag of `<density>`
 */
static std::unique_ptr<FluidProperty> createLiquidDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    config.checkConfigParameter("type", "LiquidDensity");

    //! \ogs_file_param{material__fluid__density__LiquidDensity__beta}
    const auto beta = config.getConfigParameter<double>("beta");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__rho0}
    const auto rho0 = config.getConfigParameter<double>("rho0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__temperature0}
    const auto T0 = config.getConfigParameter<double>("temperature0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__p0}
    const auto p0 = config.getConfigParameter<double>("p0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__bulk_modulus}
    const auto E = config.getConfigParameter<double>("bulk_modulus");
    return std::unique_ptr<FluidProperty>(
        new LiquidDensity(beta, rho0, T0, p0, E));
}

/**
 *     \param config  ConfigTree object which contains the input data
 *                    `<type>TemperatureDependent</type>`
 *                     and it has a tag of `<density>`
 */
static std::unique_ptr<FluidProperty> createLinearTemperatureDependentDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    config.checkConfigParameter("type", "TemperatureDependent");

    //! \ogs_file_param{material__fluid__density__TemperatureDependent__rho0}
    const auto rho0 = config.getConfigParameter<double>("rho0");
    //! \ogs_file_param{material__fluid__density__TemperatureDependent__temperature0}
    const auto T0 = config.getConfigParameter<double>("temperature0");
    //! \ogs_file_param{material__fluid__density__TemperatureDependent__beta}
    const auto beta = config.getConfigParameter<double>("beta");
    return std::unique_ptr<FluidProperty>(
        new LinearTemperatureDependentDensity(rho0, T0, beta));
}

static std::unique_ptr<FluidProperty> createLinearConcentrationDependentDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    config.checkConfigParameter("type", "ConcentrationDependent");

    const double reference_density =
    //! \ogs_file_param{material__fluid__density__ConcentrationDependent__reference_density}
        config.getConfigParameter<double>("reference_density");
    const double reference_concentration =
    //! \ogs_file_param{material__fluid__density__ConcentrationDependent__reference_concentration}
        config.getConfigParameter<double>("reference_concentration");
    const double fluid_density_difference_ratio =
    //! \ogs_file_param{material__fluid__density__ConcentrationDependent__fluid_density_difference_ratio}
        config.getConfigParameter<double>("fluid_density_difference_ratio");
    return std::unique_ptr<FluidProperty>(
        new LinearConcentrationDependentDensity(
            reference_density,
            reference_concentration,
            fluid_density_difference_ratio));
}

std::unique_ptr<FluidProperty> createFluidDensityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        //! \ogs_file_param{material__fluid__density__type}
        config.checkConfigParameter("type", "Constant");
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__density__Constant__value}
            config.getConfigParameter<double>("value")));
    }
    else if (type == "LiquidDensity")
        return createLiquidDensity(config);
    else if (type == "TemperatureDependent")
        return createLinearTemperatureDependentDensity(config);
    else if (type == "ConcentrationDependent")
        return createLinearConcentrationDependentDensity(config);
    else if (type == "IdealGasLaw")
    {
        //! \ogs_file_param{material__fluid__density__type}
        config.checkConfigParameter("type", "IdealGasLaw");
        return std::unique_ptr<FluidProperty>(
            //! \ogs_file_param{material__fluid__density__IdealGasLaw__molar_mass}
            new IdealGasLaw(config.getConfigParameter<double>("molar_mass")));
    }
    else if (type == "WaterDensityIAPWSIF97Region1")
    {
        return std::unique_ptr<FluidProperty>(
            new WaterDensityIAPWSIF97Region1());
    }
    else
    {
        OGS_FATAL(
            "The density type %s is unavailable.\n"
            "The available types are: \n\tConstant, \n\tLiquidDensity, "
            "\n\tTemperatureDependent, \n\tIdealGasLaw."
            "\n\tWaterDensityIAPWSIF97Region1\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
