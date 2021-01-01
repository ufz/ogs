/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    for (auto& phase_component : phase_components)
    {
        auto& name = phase_component.getName();
        pp_assemblage_comps[name] = *phase_component.castToBaseClass();
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
