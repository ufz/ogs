/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChemicalSolverInterface.h"
#include "Common/CreateReactionRate.h"
#include "MeshLib/Mesh.h"
#include "PhreeqcIO.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/CreateAqueousSolution.h"
#include "PhreeqcIOData/CreateEquilibriumPhase.h"
#include "PhreeqcIOData/CreateKineticReactant.h"
#include "PhreeqcIOData/CreateKnobs.h"
#include "PhreeqcIOData/CreateOutput.h"
#include "PhreeqcIOData/CreateSurface.h"
#include "PhreeqcIOData/CreateUserPunch.h"
#include "PhreeqcIOData/EquilibriumPhase.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/Knobs.h"
#include "PhreeqcIOData/ReactionRate.h"
#include "PhreeqcIOData/Surface.h"
#include "PhreeqcIOData/UserPunch.h"
#include "PhreeqcKernel.h"
#include "PhreeqcKernelData/AqueousSolution.h"
#include "PhreeqcKernelData/CreateAqueousSolution.h"
#include "PhreeqcKernelData/CreateKineticReactant.h"
#include "PhreeqcKernelData/ReactionRate.h"

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/cxxKinetics.h"

namespace
{
std::string parseDatabasePath(BaseLib::ConfigTree const& config)
{
    // database
    //! \ogs_file_param{prj__chemical_system__database}
    auto const database = config.getConfigParameter<std::string>("database");
    auto path_to_database =
        BaseLib::joinPaths(BaseLib::getProjectDirectory(), database);

    if (!BaseLib::IsFileExisting(path_to_database))
    {
        OGS_FATAL("Not found the specified thermodynamicdatabase: %s",
                  path_to_database.c_str());
    }

    INFO("Found the specified thermodynamic database: %s",
         path_to_database.c_str());

    return path_to_database;
}
}  // namespace

namespace ChemistryLib
{
template <>
std::unique_ptr<ChemicalSolverInterface>
createChemicalSolverInterface<ChemicalSolver::Phreeqc>(
    MeshLib::Mesh const& mesh,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    BaseLib::ConfigTree const& config, std::string const& output_directory)
{
    auto path_to_database = parseDatabasePath(config);

    // solution
    auto aqueous_solution = PhreeqcIOData::createAqueousSolution(
        //! \ogs_file_param{prj__chemical_system__solution}
        config.getConfigSubtree("solution"),
        process_id_to_component_name_map);

    // kinetic reactants
    auto kinetic_reactants = PhreeqcIOData::createKineticReactants(
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants}
        config.getConfigSubtreeOptional("kinetic_reactants"), mesh);

    // rates
    auto reaction_rates = createReactionRates<PhreeqcIOData::ReactionRate>(
        //! \ogs_file_param{prj__chemical_system__rates}
        config.getConfigSubtreeOptional("rates"));

    auto const num_chemical_systems = mesh.getNumberOfBaseNodes();
    std::vector<PhreeqcIOData::AqueousSolution> aqueous_solutions(
        num_chemical_systems, aqueous_solution);

    // equilibrium phases
    auto equilibrium_phases = PhreeqcIOData::createEquilibriumPhases(
        //! \ogs_file_param{prj__chemical_system__equilibrium_phases}
        config.getConfigSubtreeOptional("equilibrium_phases"), mesh);

    // surface
    auto surface = PhreeqcIOData::createSurface(
        //! \ogs_file_param{prj__chemical_system__surface}
        config.getConfigSubtreeOptional("surface"));

    // output
    auto const& components = aqueous_solution.components;
    auto const project_file_name = BaseLib::joinPaths(
        output_directory,
        BaseLib::extractBaseNameWithoutExtension(config.getProjectFileName()));

    // knobs
    auto knobs = PhreeqcIOData::createKnobs(
        //! \ogs_file_param{prj__chemical_system__knobs}
        config.getConfigSubtreeOptional("knobs"));

    // user punch
    auto user_punch = PhreeqcIOData::createUserPunch(
        //! \ogs_file_param{prj__chemical_system__user_punch}
        config.getConfigSubtreeOptional("user_punch"), mesh);
    auto output = PhreeqcIOData::createOutput(
        components, equilibrium_phases, kinetic_reactants, project_file_name);

    return std::make_unique<PhreeqcIOData::PhreeqcIO>(
        std::move(project_file_name), std::move(path_to_database),
        std::move(aqueous_solutions), std::move(equilibrium_phases),
        std::move(kinetic_reactants), std::move(reaction_rates),
        std::move(output), process_id_to_component_name_map);
}

template <>
std::unique_ptr<ChemicalSolverInterface>
createChemicalSolverInterface<ChemicalSolver::PhreeqcKernel>(
    MeshLib::Mesh const& mesh,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    BaseLib::ConfigTree const& config, std::string const& /*output_directory*/)
{
    auto path_to_database = parseDatabasePath(config);

    // solution
    auto aqueous_solution = PhreeqcKernelData::createAqueousSolution(
        //! \ogs_file_param{prj__chemical_system__solution}
        config.getConfigSubtree("solution"),
        process_id_to_component_name_map);

    // kinetic reactants
    auto kinetic_reactants = PhreeqcKernelData::createKineticReactants(
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants}
        config.getConfigSubtreeOptional("kinetic_reactants"), mesh);

    // rates
    auto reaction_rates = createReactionRates<PhreeqcKernelData::ReactionRate>(
        //! \ogs_file_param{prj__chemical_system__rates}
        config.getConfigSubtreeOptional("rates"));

    return std::make_unique<PhreeqcKernelData::PhreeqcKernel>(
        mesh.getNumberOfBaseNodes(), process_id_to_component_name_map,
        std::move(path_to_database), std::move(aqueous_solution),
        std::move(kinetic_reactants), std::move(reaction_rates));
}
}  // namespace ChemistryLib
