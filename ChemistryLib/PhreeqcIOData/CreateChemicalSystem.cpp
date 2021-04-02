/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChemicalSystem.h"

#include "BaseLib/ConfigTree.h"
#include "ChemicalSystem.h"
#include "CreateAqueousSolution.h"
#include "CreateEquilibriumReactants.h"
#include "CreateKineticReactant.h"
#include "EquilibriumReactant.h"
#include "KineticReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<ChemicalSystem> createChemicalSystem(
    BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh)
{
    // solution
    auto aqueous_solution = createAqueousSolution(
        //! \ogs_file_param{prj__chemical_system__solution}
        config.getConfigSubtree("solution"), mesh);

    // kinetic reactants
    auto kinetic_reactants = createKineticReactants(
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants}
        config.getConfigSubtreeOptional("kinetic_reactants"), mesh);

    // equilibrium reactants
    auto equilibrium_reactants = createEquilibriumReactants(
        //! \ogs_file_param{prj__chemical_system__equilibrium_reactants}
        config.getConfigSubtreeOptional("equilibrium_reactants"), mesh);

    return std::make_unique<ChemicalSystem>(std::move(aqueous_solution),
                                            std::move(kinetic_reactants),
                                            std::move(equilibrium_reactants));
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
