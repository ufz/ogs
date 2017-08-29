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
struct HTProcessData
{
    HTProcessData(
        MaterialLib::PorousMedium::PorousMediaProperties&&
            porous_media_properties_,
        ProcessLib::Parameter<double> const& density_solid_,
        ProcessLib::Parameter<double> const& fluid_reference_density_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&&
            fluid_properties_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_transversal_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_fluid_,
        Eigen::VectorXd specific_body_force_)
        : porous_media_properties(std::move(porous_media_properties_)),
          density_solid(density_solid_),
          fluid_reference_density(fluid_reference_density_),
          fluid_properties(std::move(fluid_properties_)),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          thermal_dispersivity_longitudinal(thermal_dispersivity_longitudinal_),
          thermal_dispersivity_transversal(thermal_dispersivity_transversal_),
          thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          specific_body_force(std::move(specific_body_force_))
    {
    }

    HTProcessData(HTProcessData&& other)
        : porous_media_properties(std::move(other.porous_media_properties)),
          density_solid(other.density_solid),
          fluid_reference_density(other.fluid_reference_density),
          fluid_properties(other.fluid_properties.release()),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          thermal_dispersivity_longitudinal(
              other.thermal_dispersivity_longitudinal),
          thermal_dispersivity_transversal(
              other.thermal_dispersivity_transversal),
          thermal_conductivity_solid(other.thermal_conductivity_solid),
          thermal_conductivity_fluid(other.thermal_conductivity_fluid),
          specific_body_force(other.specific_body_force)
    {
    }

    //! Copies are forbidden.
    HTProcessData(HTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HTProcessData&&) = delete;

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties;
    Parameter<double> const& density_solid;
    Parameter<double> const& fluid_reference_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_properties;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& thermal_dispersivity_longitudinal;
    Parameter<double> const& thermal_dispersivity_transversal;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
    Eigen::VectorXd const specific_body_force;
};

}  // namespace HT
}  // namespace ProcessLib
