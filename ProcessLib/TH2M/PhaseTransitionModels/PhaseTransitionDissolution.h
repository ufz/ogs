/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>

#include "MaterialLib/MPL/Medium.h"
#include "PhaseTransitionModel.h"

namespace ProcessLib
{
namespace TH2M
{
/// Dissolution only: Amounts of the gas phase can be dissolved into the liquid
/// phase. This is realized by defining two components in the liqid phase but
/// only one (or zero) components in the gas phase.
struct PhaseTransitionDissolution : PhaseTransitionModel
{
    explicit PhaseTransitionDissolution(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    PhaseTransitionModelVariables updateConstitutiveVariables(
        PhaseTransitionModelVariables const& phase_transition_model_variables,
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t,
        double const dt) const override;
};

}  // namespace TH2M
}  // namespace ProcessLib
