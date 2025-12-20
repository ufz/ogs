// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <phreeqcpp/PPassemblage.h>

#include <set>
#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class PhaseComponent final : private cxxPPassemblageComp
{
public:
    PhaseComponent(std::string&& name_, double const initial_amount,
                   double const saturation_index);

    cxxPPassemblageComp const* castToBaseClass() const
    {
        return static_cast<cxxPPassemblageComp const*>(this);
    }

    std::string const& getName() const { return Get_name(); }
};

class EquilibriumReactants final : private cxxPPassemblage
{
public:
    explicit EquilibriumReactants(
        std::vector<PhaseComponent> const& phase_components);

    void setChemicalSystemID(std::size_t const chemical_system_id)
    {
        Set_n_user_both(chemical_system_id);
    }

    cxxPPassemblage const* castToBaseClass() const
    {
        return static_cast<cxxPPassemblage const*>(this);
    }

    std::map<std::string, cxxPPassemblageComp> const& getPhaseComponents() const
    {
        return castToBaseClass()->Get_pp_assemblage_comps();
    }
};

}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
