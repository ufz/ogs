// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "InitialAqueousSolution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
InitialAqueousSolution::InitialAqueousSolution(
    std::map<std::string, cxxISolutionComp>& components)
{
    units = "Mol/kgw";

    for (auto& map_pair : components)
    {
        auto& component = map_pair.second;
        component.Set_units(units.c_str());
        component.Set_pe_reaction(default_pe);
    }

    comps = components;
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
