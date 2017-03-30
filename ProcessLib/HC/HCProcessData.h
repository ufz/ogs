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

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "PorousMediaProperties.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace HC
{
struct HCProcessData
{
    HCProcessData(
        PorousMediaProperties&& porous_media_properties_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity_model_,
        ProcessLib::Parameter<double> const& fluid_reference_density_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& fluid_density_,
        ProcessLib::Parameter<double> const& molecular_diffusion_coefficient_,
        ProcessLib::Parameter<double> const& solute_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& solute_dispersivity_transverse_,
        ProcessLib::Parameter<double> const& retardation_factor_,
        ProcessLib::Parameter<double> const& decay_rate_,
        Eigen::VectorXd const& specific_body_force_,
        bool const has_gravity_)
        : porous_media_properties(std::move(porous_media_properties_)),
          viscosity_model(std::move(viscosity_model_)),
          fluid_reference_density(fluid_reference_density_),
          fluid_density(std::move(fluid_density_)),
          molecular_diffusion_coefficient(molecular_diffusion_coefficient_),
          solute_dispersivity_longitudinal(solute_dispersivity_longitudinal_),
          solute_dispersivity_transverse(solute_dispersivity_transverse_),
          retardation_factor(retardation_factor_),
          decay_rate(decay_rate_),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_)
    {
    }

    HCProcessData(HCProcessData&& other)
        : porous_media_properties(std::move(other.porous_media_properties)),
          viscosity_model(other.viscosity_model.release()),
          fluid_reference_density(other.fluid_reference_density),
          fluid_density(other.fluid_density.release()),
          molecular_diffusion_coefficient(
              other.molecular_diffusion_coefficient),
          solute_dispersivity_longitudinal(
              other.solute_dispersivity_longitudinal),
          solute_dispersivity_transverse(other.solute_dispersivity_transverse),
          retardation_factor(other.retardation_factor),
          decay_rate(other.decay_rate),
          specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity)
    {
    }

    //! Copies are forbidden.
    HCProcessData(HCProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HCProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HCProcessData&&) = delete;

    PorousMediaProperties porous_media_properties;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> viscosity_model;
    Parameter<double> const& fluid_reference_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> fluid_density;
    Parameter<double> const& molecular_diffusion_coefficient;
    Parameter<double> const& solute_dispersivity_longitudinal;
    Parameter<double> const& solute_dispersivity_transverse;
    Parameter<double> const& retardation_factor;
    Parameter<double> const& decay_rate;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace HC
}  // namespace ProcessLib
