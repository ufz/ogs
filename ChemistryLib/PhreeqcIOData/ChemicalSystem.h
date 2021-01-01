/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "AqueousSolution.h"
#include "EquilibriumReactant.h"
#include "KineticReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ChemicalSystem
{
    ChemicalSystem(std::unique_ptr<AqueousSolution>&& aqueous_solution_,
                   std::vector<KineticReactant>&& kinetic_reactants_,
                   std::vector<EquilibriumReactant>&& equilibrium_reactants_)
        : aqueous_solution(std::move(aqueous_solution_)),
          kinetic_reactants(std::move(kinetic_reactants_)),
          equilibrium_reactants(std::move(equilibrium_reactants_))
    {
    }

    void initialize(std::size_t const num_chemical_systems);

    std::unique_ptr<AqueousSolution> aqueous_solution;
    std::vector<KineticReactant> kinetic_reactants;
    std::vector<EquilibriumReactant> equilibrium_reactants;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
