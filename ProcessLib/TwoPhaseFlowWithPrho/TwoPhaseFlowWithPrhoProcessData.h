/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    Eigen::VectorXd const specific_body_force_;

    bool const has_gravity_;
    bool const has_mass_lumping_;
    ParameterLib::Parameter<double> const& diffusion_coeff_component_b_;
    ParameterLib::Parameter<double> const& diffusion_coeff_component_a_;
    ParameterLib::Parameter<double> const& temperature_;
    std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties> material_;
};

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
