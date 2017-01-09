/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace HT
{
struct HTProcessData
{
    HTProcessData(
        ProcessLib::Parameter<double> const& porosity_,
        ProcessLib::Parameter<double> const& intrinsic_permeability_,
        ProcessLib::Parameter<double> const& specific_storage_,
        ProcessLib::Parameter<double> const& viscosity_,
        ProcessLib::Parameter<double> const& density_solid_,
        ProcessLib::Parameter<double> const& density_fluid_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& thermal_dispersivity_transversal_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_solid_,
        ProcessLib::Parameter<double> const& specific_heat_capacity_fluid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_solid_,
        ProcessLib::Parameter<double> const& thermal_conductivity_fluid_,
        ProcessLib::Parameter<double> const& thermal_expansion_coefficient_,
        ProcessLib::Parameter<double> const&
            reference_temperature_fluid_density_model_,
        Eigen::Vector3d const& specific_body_force_,
        bool const has_gravity_)
        : porosity(porosity_),
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          viscosity(viscosity_),
          density_solid(density_solid_),
          density_fluid(density_fluid_),
          specific_heat_capacity_solid(specific_heat_capacity_solid_),
          specific_heat_capacity_fluid(specific_heat_capacity_fluid_),
          thermal_dispersivity_longitudinal(thermal_dispersivity_longitudinal_),
          thermal_dispersivity_transversal(thermal_dispersivity_transversal_),
          thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          thermal_expansion_coefficient(thermal_expansion_coefficient_),
          reference_temperature_fluid_density_model(
              reference_temperature_fluid_density_model_),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_)
    {
    }

    HTProcessData(HTProcessData&& other)
        : porosity(other.porosity),
          intrinsic_permeability(other.intrinsic_permeability),
          specific_storage(other.specific_storage),
          viscosity(other.viscosity),
          density_solid(other.density_solid),
          density_fluid(other.density_fluid),
          specific_heat_capacity_solid(other.specific_heat_capacity_solid),
          specific_heat_capacity_fluid(other.specific_heat_capacity_fluid),
          thermal_dispersivity_longitudinal(other.thermal_dispersivity_longitudinal),
          thermal_dispersivity_transversal(other.thermal_dispersivity_transversal),
          thermal_conductivity_solid(other.thermal_conductivity_solid),
          thermal_conductivity_fluid(other.thermal_conductivity_fluid),
          thermal_expansion_coefficient(other.thermal_expansion_coefficient),
          reference_temperature_fluid_density_model(
              other.reference_temperature_fluid_density_model),
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

    Parameter<double> const& porosity;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& specific_storage;
    Parameter<double> const& viscosity;
    Parameter<double> const& density_solid;
    Parameter<double> const& density_fluid;
    Parameter<double> const& specific_heat_capacity_solid;
    Parameter<double> const& specific_heat_capacity_fluid;
    Parameter<double> const& thermal_dispersivity_longitudinal;
    Parameter<double> const& thermal_dispersivity_transversal;
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
    Parameter<double> const& thermal_expansion_coefficient;
    Parameter<double> const& reference_temperature_fluid_density_model;
    Eigen::Vector3d const specific_body_force;
    bool const has_gravity;
};

}  // namespace HT
}  // namespace ProcessLib
