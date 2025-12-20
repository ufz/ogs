// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EquilibriumReactants.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
PhaseComponent::PhaseComponent(std::string&& name_,
                               double const initial_amount,
                               double const saturation_index)
{
    name = std::move(name_);
    moles = initial_amount;
    si = saturation_index;
    si_org = saturation_index;
}

EquilibriumReactants::EquilibriumReactants(
    std::vector<PhaseComponent> const& phase_components)
{
    for (auto const& phase_component : phase_components)
    {
        auto& name = phase_component.getName();
        pp_assemblage_comps[name] = *phase_component.castToBaseClass();
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
