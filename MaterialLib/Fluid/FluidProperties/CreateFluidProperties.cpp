/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateFluidProperties.cpp
 *
 * Created on December 13, 2016, 3:32 PM
 */

#include "CreateFluidProperties.h"

#include <string>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"

#include "FluidProperties.h"
#include "PrimaryVariableDependentFluidProperties.h"
#include "FluidPropertiesWithDensityDependentModels.h"

namespace MaterialLib
{
namespace Fluid
{
std::unique_ptr<FluidProperties> createFluidProperties(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fluid__density}
    auto const& rho_conf = config.getConfigSubtree("density");
    auto liquid_density = MaterialLib::Fluid::createFluidDensityModel(rho_conf);

    //! \ogs_file_param{material__fluid__viscosity}
    auto const& mu_conf = config.getConfigSubtree("viscosity");
    auto viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);
    const bool is_mu_density_dependent =
        (viscosity->getName().find("density dependent") != std::string::npos);

    bool is_cp_density_dependent = false;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> specific_heat_capacity =
        nullptr;
    auto heat_capacity__opt_conf =
        //! \ogs_file_param{material__fluid__specific_heat_capacity}
        config.getConfigSubtreeOptional("specific_heat_capacity");
    if (heat_capacity__opt_conf)
    {
        const auto& heat_capacity_conf = *heat_capacity__opt_conf;
        specific_heat_capacity =
            createSpecificFluidHeatCapacityModel(heat_capacity_conf);
        is_cp_density_dependent =
            (specific_heat_capacity->getName().find("density dependent") !=
             std::string::npos);
    }

    bool is_KT_density_dependent = false;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> thermal_conductivity =
        nullptr;
    auto const& thermal_conductivity_opt_conf =
        //! \ogs_file_param{material__fluid__thermal_conductivity}
        config.getConfigSubtreeOptional("thermal_conductivity");
    if (thermal_conductivity_opt_conf)
    {
        auto const& thermal_conductivity_conf = *thermal_conductivity_opt_conf;
        thermal_conductivity =
            MaterialLib::Fluid::createFluidThermalConductivityModel(
                thermal_conductivity_conf);
        is_KT_density_dependent =
            (specific_heat_capacity->getName().find("density dependent") !=
             std::string::npos);
    }

    if (is_mu_density_dependent || is_cp_density_dependent ||
        is_KT_density_dependent)
        return std::make_unique<
            MaterialLib::Fluid::FluidPropertiesWithDensityDependentModels>(
            std::move(liquid_density), std::move(viscosity),
            std::move(specific_heat_capacity), std::move(thermal_conductivity),
            is_mu_density_dependent, is_cp_density_dependent,
            is_KT_density_dependent);

    return std::make_unique<
        MaterialLib::Fluid::PrimaryVariableDependentFluidProperties>(
        std::move(liquid_density), std::move(viscosity),
        std::move(specific_heat_capacity), std::move(thermal_conductivity));
}

}  // end namespace
}  // end namespace
