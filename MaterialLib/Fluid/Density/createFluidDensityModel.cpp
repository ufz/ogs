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

#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "IdealGasLaw.h"
#include "LinearTemperatureDependentDensity.h"
#include "LiquidDensity.h"

namespace MaterialLib
{
namespace Fluid
{
/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>fluid</type> and it has
                   a tag of <density>
*/
static std::unique_ptr<FluidProperty> createLiquidDensity(BaseLib::ConfigTree const& config)
{
    std::array<double, 5> parameters = {
        {//! \ogs_file_param{material__fluid__density__liquid_density__beta}
         config.getConfigParameter<double>("beta"),
         //! \ogs_file_param{material__fluid__density__liquid_density__rho0}
         config.getConfigParameter<double>("rho0"),
         //! \ogs_file_param{material__fluid__density__liquid_density__temperature0}
         config.getConfigParameter<double>("temperature0"),
         //! \ogs_file_param{material__fluid__density__liquid_density__p0}
         config.getConfigParameter<double>("p0"),
         //! \ogs_file_param{material__fluid__density__liquid_density__bulk_modulus}
         config.getConfigParameter<double>("bulk_modulus")}};
    return std::unique_ptr<FluidProperty>(new LiquidDensity(parameters));
}

/*!
    \param config  ConfigTree object which contains the input data
                   including  <type>fluid</type> and it has
                   a tag of <density>
*/
static std::unique_ptr<FluidProperty> createLinearTemperatureDependentDensity
                                     (BaseLib::ConfigTree const& config)
{
    std::array<double, 3> parameters = {
        {//! \ogs_file_param{material__fluid__density__linear_temperature__rho0}
         config.getConfigParameter<double>("rho0"),
         //! \ogs_file_param{material__fluid__density__linear_temperature__temperature0}
         config.getConfigParameter<double>("temperature0"),
         //! \ogs_file_param{material__fluid__density__linear_temperature__beta}
         config.getConfigParameter<double>("beta")}};
    return std::unique_ptr<FluidProperty>(
            new LinearTemperatureDependentDensity(parameters));
}

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
        return createLiquidDensity(config);

    //! \ogs_file_param{material__fluid__density__TemperatureDependent}
    else if (type == "TemperatureDependent")
        return createLinearTemperatureDependentDensity(config);
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
            "The density type %s is unavailable.\n"
            "The available types are: \n\tConstant, \n\tLiquidDensity, "
            "\n\tTemperatureDependent, \n\tIdealGasLaw.\n", type.data());
    }
}

}  // end namespace
}  // end namespace
