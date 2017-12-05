/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/PorousMedium/PorousMediaProperties.h"

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace HT
{
struct HTMaterialProperties
{
    HTMaterialProperties(
        MaterialLib::PorousMedium::PorousMediaProperties&&
            porous_media_properties_,
        ProcessLib::Parameter<double> const& density_solid_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&&
            fluid_properties_,
        bool const has_fluid_thermal_dispersivity_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_transversal_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_fluid_,
        bool const has_fluid_thermal_expansion_,
        Parameter<double> const& solid_thermal_expansion_,
        Parameter<double> const& biot_constant_,
        Eigen::VectorXd specific_body_force_,
        bool const has_gravity_)
        : porous_media_properties(std::move(porous_media_properties_)),
          density_solid(density_solid_),
          fluid_properties(std::move(fluid_properties_)),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          has_fluid_thermal_dispersivity(has_fluid_thermal_dispersivity_),
          thermal_dispersivity_longitudinal(thermal_dispersivity_longitudinal_),
          thermal_dispersivity_transversal(thermal_dispersivity_transversal_),
          thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          has_fluid_thermal_expansion(has_fluid_thermal_expansion_),
          solid_thermal_expansion(solid_thermal_expansion_),
          biot_constant(biot_constant_),
          specific_body_force(std::move(specific_body_force_)),
          has_gravity(has_gravity_)
    {
    }

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties;
    Parameter<double> const& density_solid;
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_properties;
    Parameter<double> const& specific_heat_capacity_solid;
    bool const has_fluid_thermal_dispersivity;
    Parameter<double> const& thermal_dispersivity_longitudinal;
    Parameter<double> const& thermal_dispersivity_transversal;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;

    bool const has_fluid_thermal_expansion;
    Parameter<double> const& solid_thermal_expansion;
    Parameter<double> const& biot_constant;

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace HT
}  // namespace ProcessLib
