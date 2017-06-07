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

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "PorousMediaProperties.h"

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace HT
{
struct HTProcessData
{
    HTProcessData(
        PorousMediaProperties&& porous_media_properties_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity_model_,
        ProcessLib::Parameter<double> const& density_solid_,
        ProcessLib::Parameter<double> const& fluid_reference_density_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& fluid_density_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_transversal_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_fluid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_fluid_,
        Eigen::VectorXd specific_body_force_,
        bool const has_gravity_)
        : porous_media_properties(std::move(porous_media_properties_)),
          viscosity_model(std::move(viscosity_model_)),
          density_solid(density_solid_),
          fluid_reference_density(fluid_reference_density_),
          fluid_density(std::move(fluid_density_)),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          specific_heat_capacity_fluid(specific_heat_capacity_fluid_),
          thermal_dispersivity_longitudinal(thermal_dispersivity_longitudinal_),
          thermal_dispersivity_transversal(thermal_dispersivity_transversal_),
          thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          specific_body_force(std::move(specific_body_force_)),
          has_gravity(has_gravity_)
    {
    }

    HTProcessData(HTProcessData&& other)
        : porous_media_properties(std::move(other.porous_media_properties)),
          viscosity_model(other.viscosity_model.release()),
          density_solid(other.density_solid),
          fluid_reference_density(other.fluid_reference_density),
          fluid_density(other.fluid_density.release()),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          specific_heat_capacity_fluid(other.specific_heat_capacity_fluid),
          thermal_dispersivity_longitudinal(
              other.thermal_dispersivity_longitudinal),
          thermal_dispersivity_transversal(
              other.thermal_dispersivity_transversal),
          thermal_conductivity_solid(other.thermal_conductivity_solid),
          thermal_conductivity_fluid(other.thermal_conductivity_fluid),
          specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity)
    {
    }

    //! Copies are forbidden.
    HTProcessData(HTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HTProcessData&&) = delete;

    PorousMediaProperties porous_media_properties;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> viscosity_model;
    Parameter<double> const& density_solid;
    Parameter<double> const& fluid_reference_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> fluid_density;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& specific_heat_capacity_fluid;
    Parameter<double> const& thermal_dispersivity_longitudinal;
    Parameter<double> const& thermal_dispersivity_transversal;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace HT
}  // namespace ProcessLib
