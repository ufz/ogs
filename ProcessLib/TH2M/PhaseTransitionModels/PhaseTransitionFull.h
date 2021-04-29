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
#include "PhaseTransitionModels.h"

namespace ProcessLib
{
namespace TH2M
{
/// Full phase transition: Gas can dissolve into the liquid phase according to
/// an equilibrium and water evaporates into the gas phase according to Dalton's
/// law. This is realized by defining two components in each gas and liquid
/// phase.
struct PhaseTransitionFull : PhaseTransitionModels
{
    explicit PhaseTransitionFull(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void computeConstitutiveVariables(
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t,
        double const dt) override;
};

}  // namespace TH2M
}  // namespace ProcessLib
