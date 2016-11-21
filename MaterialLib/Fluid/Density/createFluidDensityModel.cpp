/*!
   \file  createFluidDensityModel.cpp
   \brief create an instance of a fluid density class.

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <array>

#include "createFluidDensityModel.h"

#include "BaseLib/Error.h"

#include "IdealGasLaw.h"
#include "LinearTemperatureDependentDensity.h"
#include "LiquidDensity.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>LiquidDensity</type> and it has
                   a tag of <density>
*/
static std::unique_ptr<FluidProperty> createLiquidDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__LiquidDensity}
    config.checkConfigParameter("type", "LiquidDensity");

    //! \ogs_file_param{material__fluid__density__LiquidDensity__beta}
    const double beta = config.getConfigParameter<double>("beta");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__rho0}
    const double rho0 = config.getConfigParameter<double>("rho0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__temperature0}
    const double T0 = config.getConfigParameter<double>("temperature0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__p0}
    const double p0 = config.getConfigParameter<double>("p0");
    //! \ogs_file_param{material__fluid__density__LiquidDensity__bulk_modulus}
    const double E = config.getConfigParameter<double>("bulk_modulus");
    return std::unique_ptr<FluidProperty>(
        new LiquidDensity(beta, rho0, T0, p0, E));
}

/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>TemperatureDependent</type> and it has
                   a tag of <density>
*/
static std::unique_ptr<FluidProperty> createLinearTemperatureDependentDensity(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__TemperatureDependent}
    config.checkConfigParameter("type", "TemperatureDependent");

    //! \ogs_file_param{material__fluid__density__TemperatureDependent__rho0}
    const double rho0 = config.getConfigParameter<double>("rho0");
    //! \ogs_file_param{material__fluid__density__TemperatureDependent__temperature0}
    const double T0 = config.getConfigParameter<double>("temperature0");
    //! \ogs_file_param{material__fluid__density__TemperatureDependent__beta}
    const double beta = config.getConfigParameter<double>("beta");
    return std::unique_ptr<FluidProperty>(
        new LinearTemperatureDependentDensity(rho0, T0, beta));
}

std::unique_ptr<FluidProperty> createFluidDensityModel(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        //! \ogs_file_param{material__fluid__density__Constant}
        config.checkConfigParameter("type", "Constant");
        return std::unique_ptr<FluidProperty>(new ConstantFluidProperty(
            //! \ogs_file_param{material__fluid__density__Constant__value}
            config.getConfigParameter<double>("value")));
    }
    else if (type == "LiquidDensity")
        return createLiquidDensity(config);
    else if (type == "TemperatureDependent")
        return createLinearTemperatureDependentDensity(config);
    else if (type == "IdealGasLaw")
    {
        //! \ogs_file_param{material__fluid__density__IdealGasLaw}
        config.checkConfigParameter("type", "IdealGasLaw");
        return std::unique_ptr<FluidProperty>(
            //! \ogs_file_param{material__fluid__density__IdealGasLaw__molar_mass}
            new IdealGasLaw(config.getConfigParameter<double>("molar_mass")));
    }
    else
    {
        OGS_FATAL(
            "The density type %s is unavailable.\n"
            "The available types are: \n\tConstant, \n\tLiquidDensity, "
            "\n\tTemperatureDependent, \n\tIdealGasLaw.\n",
            type.data());
    }
}

}  // end namespace
}  // end namespace
