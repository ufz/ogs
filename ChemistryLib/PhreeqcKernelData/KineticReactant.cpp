// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "KineticReactant.h"

#include <algorithm>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
KineticReactant::KineticReactant(std::string name, double const initial_amount)
{
    rate_name = std::move(name);
    namecoef.add(rate_name.c_str(), 1.0);
    m = initial_amount;
    m0 = initial_amount;
}

Kinetics::Kinetics(std::vector<KineticReactant> const& kinetic_reactants)
{
    std::transform(kinetic_reactants.begin(),
                   kinetic_reactants.end(),
                   std::back_inserter(kinetics_comps),
                   [](KineticReactant const& kinetic_reactant)
                   { return *kinetic_reactant.castToBaseClass(); });
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
