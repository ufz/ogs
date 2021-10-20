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
#include <variant>
#include <vector>

#include "AqueousSolution.h"
#include "EquilibriumReactant.h"
#include "Exchange.h"
#include "KineticReactant.h"
#include "Surface.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ChemicalSystem
{
    ChemicalSystem(std::unique_ptr<AqueousSolution>&& aqueous_solution_,
                   std::vector<KineticReactant>&& kinetic_reactants_,
                   std::vector<EquilibriumReactant>&& equilibrium_reactants_,
                   std::vector<ExchangeSite>&& exchangers_,
                   std::vector<std::variant<DensityBasedSurfaceSite,
                                            MoleBasedSurfaceSite>>&& surface_)
        : aqueous_solution(std::move(aqueous_solution_)),
          kinetic_reactants(std::move(kinetic_reactants_)),
          equilibrium_reactants(std::move(equilibrium_reactants_)),
          exchangers(std::move(exchangers_)),
          surface(std::move(surface_))
    {
    }

    void initialize(std::size_t const num_chemical_systems);

    std::unique_ptr<AqueousSolution> aqueous_solution;
    std::vector<KineticReactant> kinetic_reactants;
    std::vector<EquilibriumReactant> equilibrium_reactants;
    std::vector<ExchangeSite> exchangers;
    std::vector<std::variant<DensityBasedSurfaceSite, MoleBasedSurfaceSite>>
        surface;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
