/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ChemicalSystem.h"

#include "AqueousSolution.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void ChemicalSystem::initialize(std::size_t const num_chemical_systems)
{
    aqueous_solution->pH =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            num_chemical_systems);

    aqueous_solution->pe->resize(num_chemical_systems, aqueous_solution->pe0);

    auto& components = aqueous_solution->components;
    for (auto& component : components)
    {
        component.amount =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                num_chemical_systems);
    }

    for (auto& kinetic_reactant : kinetic_reactants)
    {
        kinetic_reactant.molality->resize(num_chemical_systems);
        kinetic_reactant.molality_prev->resize(num_chemical_systems);
        kinetic_reactant.volume_fraction->resize(num_chemical_systems);
        kinetic_reactant.volume_fraction_prev->resize(num_chemical_systems);
    }

    for (auto& equilibrium_reactant : equilibrium_reactants)
    {
        equilibrium_reactant.molality->resize(num_chemical_systems);
        equilibrium_reactant.molality_prev->resize(num_chemical_systems);
        equilibrium_reactant.volume_fraction->resize(num_chemical_systems);
        equilibrium_reactant.volume_fraction_prev->resize(num_chemical_systems);
    }

    for (auto& exchanger : exchangers)
    {
        exchanger.molality->resize(num_chemical_systems);
    }

    for (auto& surface_site : surface)
    {
        if (auto const surface_site_ptr =
                std::get_if<MoleBasedSurfaceSite>(&surface_site))
        {
            surface_site_ptr->molality->resize(num_chemical_systems);
        }
    }
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
