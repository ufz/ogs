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
struct PhaseTransitionEvaporation : PhaseTransitionModel
{
    explicit PhaseTransitionEvaporation(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void computeConstitutiveVariables(
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t,
        double const dt) override;

private:
    int const n_components_gas_;
    int const gas_phase_vapour_component_index_;
    int const gas_phase_dry_air_component_index_;
};

}  // namespace TH2M
}  // namespace ProcessLib
