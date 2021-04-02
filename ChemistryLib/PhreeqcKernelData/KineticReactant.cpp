/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "KineticReactant.h"

#include <algorithm>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
KineticReactant::KineticReactant(std::string const& name,
                                 double const initial_amount)
{
    rate_name = name.c_str();
    namecoef.add(name.c_str(), 1.0);
    m = initial_amount;
    m0 = initial_amount;
}

Kinetics::Kinetics(std::vector<KineticReactant> const& kinetic_reactants)
{
    std::transform(kinetic_reactants.begin(),
                   kinetic_reactants.end(),
                   std::back_inserter(kinetics_comps),
                   [](KineticReactant const& kinetic_reactant) {
                       return *kinetic_reactant.castToBaseClass();
                   });
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
