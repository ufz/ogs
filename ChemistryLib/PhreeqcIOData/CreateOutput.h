/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "AqueousSolution.h"
#include "EquilibriumPhase.h"
#include "KineticReactant.h"
#include "Output.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<Output> createOutput(
    std::vector<Component> const& components,
    std::vector<EquilibriumPhase> const& equilibrium_phases,
    std::vector<KineticReactant> const& kinetic_reactants,
    std::string const& project_file_name);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
