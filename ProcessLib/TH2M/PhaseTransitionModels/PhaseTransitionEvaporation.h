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
struct PhaseTransitionEvaporation : PhaseTransitionModels
{
    explicit PhaseTransitionEvaporation(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media);

    void getConstitutiveVariables(const MaterialPropertyLib::Medium* medium,
                                  MaterialPropertyLib::VariableArray variables,
                                  ParameterLib::SpatialPosition pos,
                                  double const t, double const dt) override;

private:
    int n_components_gas_;
    int gas_phase_vapour_component_index_ = -1;
    int gas_phase_dry_air_component_index_ = -1;
};

}  // namespace TH2M
}  // namespace ProcessLib
