// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <variant>
#include <vector>

#include "AqueousSolution.h"
#include "EquilibriumReactant.h"
#include "Exchange.h"
#include "KineticReactant.h"
#include "Surface.h"

/**
 * \file
 * \brief Definition of one reactive chemical system for PHREEQC coupling.
 *
 * OpenGeoSys treats each reactive control volume as an independent chemical
 * system when calling PHREEQC. That local system includes not only the
 * aqueous fluid, but also solid and sorbed phases associated with that
 * volume (kinetic reactants, equilibrium minerals, ion exchangers, and
 * surface complexation sites).
 *
 * For each \c chemical_system_id, PHREEQC is provided with:
 *  - the aqueous state (AqueousSolution: component totals, pH / pe control,
 *    temperature, pressure),
 *  - kinetic reactants (rate-controlled mass transfer),
 *  - equilibrium reactants (phases enforced at equilibrium),
 *  - ion-exchange sites (ExchangeSite),
 *  - surface complexation / sorption sites (DensityBasedSurfaceSite,
 *    MoleBasedSurfaceSite).
 *
 * PHREEQC advances this local system as a closed, well-mixed batch reactor
 * over the current time increment. After reaction, updated aqueous totals,
 * pH, mineral amounts, etc. are written back into OpenGeoSys and become the
 * starting state for the next transport step.
 *
 * initialize(num_chemical_systems) prepares internal storage so that
 * per-system quantities (e.g. per-component totals, pH) are allocated for
 * all chemical_system_id in the mesh.
 */
namespace ChemistryLib
{
namespace PhreeqcIOData
{
/**
 * \struct ChemicalSystem
 * \brief Complete description of one local reactive system passed to PHREEQC.
 *
 * For each \c chemical_system_id, ChemicalSystem groups the aqueous solution
 * state, kinetic reactants, equilibrium reactants, ion exchangers, and surface
 * complexation sites. PhreeqcIO uses this object to build PHREEQC input for
 * each \c chemical_system_id, run PHREEQC, and collect the reacted state.
 *
 * initialize(num_chemical_systems) resizes/initializes internal data so that
 * all per-system vectors are ready for every \c chemical_system_id.
 */
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
