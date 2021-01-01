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

#include "TwoPhaseFlowWithPrhoMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPrho
{
struct TwoPhaseFlowWithPrhoProcessData
{
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;
    bool const _has_mass_lumping;
    ParameterLib::Parameter<double> const& _diffusion_coeff_component_b;
    ParameterLib::Parameter<double> const& _diffusion_coeff_component_a;
    ParameterLib::Parameter<double> const& _temperature;
    std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties> _material;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
