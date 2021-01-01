/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermalTwoPhaseFlowWithPP
{
struct ThermalTwoPhaseFlowWithPPProcessData
{
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;
    bool const has_mass_lumping;
    ParameterLib::Parameter<double> const& diffusion_coeff_component_b;
    ParameterLib::Parameter<double> const& diffusion_coeff_component_a;
    ParameterLib::Parameter<double> const& density_solid;
    ParameterLib::Parameter<double> const& latent_heat_evaporation;
    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties> material;
};

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
