/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/Error.h"
#include "CreateOutput.h"
#include "CreatePhreeqcIO.h"
#include "PhreeqcIO.h"
#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/CreateAqueousSolution.h"
#include "PhreeqcIOData/CreateEquilibriumPhase.h"
#include "PhreeqcIOData/CreateKineticReactant.h"
#include "PhreeqcIOData/CreateReactionRate.h"
#include "PhreeqcIOData/EquilibriumPhase.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/ReactionRate.h"

namespace ChemistryLib
{
std::unique_ptr<PhreeqcIO> createPhreeqcIO(
    std::size_t const num_nodes,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    BaseLib::ConfigTree const& config,
    std::string const& output_directory)
{
    auto const& num_chemical_systems = num_nodes;

    // database
    //! \ogs_file_param{prj__chemical_system__database}
    auto const database = config.getConfigParameter<std::string>("database");
    auto path_to_database =
        BaseLib::joinPaths(BaseLib::getProjectDirectory(), database);

    if (BaseLib::IsFileExisting(path_to_database))
    {
        INFO("Found the specified thermodynamic database: %s",
             path_to_database.c_str());
    }
    else
    {
        OGS_FATAL("Not found the specified thermodynamicdatabase: %s",
                  path_to_database.c_str());
    }

    // solution
    auto const aqueous_solution_per_chem_sys = createAqueousSolution(
        //! \ogs_file_param{prj__chemical_system__solution}
        config.getConfigSubtree("solution"));

    auto const& components_per_chem_sys =
        aqueous_solution_per_chem_sys.components;
    for (auto const& component : components_per_chem_sys)
    {
        auto process_id_to_component_name_map_element = std::find_if(
            process_id_to_component_name_map.begin(),
            process_id_to_component_name_map.end(),
            [&component](std::pair<int, std::string> const& map_element) {
                return map_element.second == component.name;
            });

        if (process_id_to_component_name_map_element ==
            process_id_to_component_name_map.end())
        {
            OGS_FATAL(
                "Component %s given in <solution>/<components> is not found in "
                "specified coupled processes (see "
                "<process>/<process_variables>/<concentration>).",
                component.name.c_str());
        }
    }
    if (components_per_chem_sys.size() + 1 !=
        process_id_to_component_name_map.size())
    {
        OGS_FATAL(
            "The number of components given in <solution>/<components> is not "
            "in line with the number of transport processes - 1 which stands "
            "for the transport process of hydrogen.");
    }

    std::vector<AqueousSolution> aqueous_solutions(
        num_chemical_systems, aqueous_solution_per_chem_sys);

    // equilibrium phases
    auto const equilibrium_phases_per_chem_sys = createEquilibriumPhases(
        //! \ogs_file_param{prj__chemical_system__equilibrium_phases}
        config.getConfigSubtreeOptional("equilibrium_phases"));
    std::vector<std::vector<EquilibriumPhase>> equilibrium_phases(
        num_chemical_systems, equilibrium_phases_per_chem_sys);

    // kinetic reactants
    auto const kinetic_reactants_per_chem_sys = createKineticReactants(
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants}
        config.getConfigSubtreeOptional("kinetic_reactants"));
    std::vector<std::vector<KineticReactant>> kinetic_reactants(
        num_chemical_systems, kinetic_reactants_per_chem_sys);

    // rates
    auto reaction_rates = createReactionRates(
        //! \ogs_file_param{prj__chemical_system__rates}
        config.getConfigSubtreeOptional("rates"));

    // output
    auto const project_file_name = BaseLib::joinPaths(
        output_directory,
        BaseLib::extractBaseNameWithoutExtension(config.getProjectFileName()));
    auto output =
        createOutput(components_per_chem_sys, equilibrium_phases_per_chem_sys,
                     kinetic_reactants_per_chem_sys, project_file_name);

    return std::make_unique<PhreeqcIO>(
        project_file_name, std::move(path_to_database),
        std::move(aqueous_solutions), std::move(equilibrium_phases),
        std::move(kinetic_reactants), std::move(reaction_rates),
        std::move(output), process_id_to_component_name_map);
}
}  // namespace ChemistryLib
