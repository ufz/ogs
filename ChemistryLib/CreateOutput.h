/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "Output.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/EquilibriumPhase.h"
#include "PhreeqcIOData/KineticReactant.h"

namespace ChemistryLib
{
std::unique_ptr<Output> createOutput(
    std::vector<Component> const& components,
    std::vector<EquilibriumPhase> const& equilibrium_phases,
    std::vector<KineticReactant> const& kinetic_reactants,
    std::string const& project_file_name);
}  // namespace ChemistryLib
