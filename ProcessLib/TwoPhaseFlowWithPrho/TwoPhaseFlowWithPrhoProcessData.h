/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TwoPhaseFlowWithPrhoMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPrho
{
struct TwoPhaseFlowWithPrhoProcessData
{
    TwoPhaseFlowWithPrhoProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& diffusion_coeff_component_b_,
        Parameter<double> const& diffusion_coeff_component_a_,
        Parameter<double> const& temperature_,
        std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>&& material_)
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
          _diffusion_coeff_component_b(diffusion_coeff_component_b_),
          _diffusion_coeff_component_a(diffusion_coeff_component_a_),
          _temperature(temperature_),
          _material(std::move(material_))

    {
    }

    TwoPhaseFlowWithPrhoProcessData(TwoPhaseFlowWithPrhoProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
          _diffusion_coeff_component_b(other._diffusion_coeff_component_b),
          _diffusion_coeff_component_a(other._diffusion_coeff_component_a),
          _temperature(other._temperature),
          _material(std::move(other._material))
    {
    }

    //! Copies are forbidden.
    TwoPhaseFlowWithPrhoProcessData(TwoPhaseFlowWithPrhoProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPrhoProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPrhoProcessData&&) = delete;
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;
    bool const _has_mass_lumping;
    Parameter<double> const& _diffusion_coeff_component_b;
    Parameter<double> const& _diffusion_coeff_component_a;
    Parameter<double> const& _temperature;
    std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties> _material;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
