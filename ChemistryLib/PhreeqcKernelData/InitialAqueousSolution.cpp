/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
