/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <cstdlib>

#include "BaseLib/Error.h"
#include "PhreeqcKernel.h"
#include "PhreeqcKernelData/AqueousSolution.h"
#include "PhreeqcKernelData/ReactionRate.h"

#include <iphreeqc/src/src/phreeqcpp/cxxKinetics.h>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
PhreeqcKernel::PhreeqcKernel(
    std::size_t const num_chemical_systems,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    std::string const& database,
    AqueousSolution aqueous_solution,
    std::unique_ptr<EquilibriumReactants>&& equilibrium_reactants,
    std::unique_ptr<Kinetics>&& kinetic_reactants,
    std::vector<ReactionRate>&& reaction_rates)
    : _initial_aqueous_solution(aqueous_solution.getInitialAqueousSolution()),
      _aqueous_solution(aqueous_solution.castToBaseClassNoninitialized()),
      _reaction_rates(std::move(reaction_rates))
{
    initializePhreeqcGeneralSettings();

    loadDatabase(database);

    // solution
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        aqueous_solution.setChemicalSystemID(chemical_system_id);
        Rxn_solution_map.emplace(chemical_system_id,
                                 *aqueous_solution.castToBaseClass());
    }
    use.Set_solution_in(true);

    // equilibrium reactants
    if (equilibrium_reactants)
    {
        tidyEquilibriumReactants(*equilibrium_reactants);

        for (std::size_t chemical_system_id = 0;
             chemical_system_id < num_chemical_systems;
             ++chemical_system_id)
        {
            equilibrium_reactants->setChemicalSystemID(chemical_system_id);
            Rxn_pp_assemblage_map.emplace(
                chemical_system_id, *equilibrium_reactants->castToBaseClass());
        }
        // explicitly release and delete equilibrium_reactants
        equilibrium_reactants.reset(nullptr);

        use.Set_pp_assemblage_in(true);
    }

    // kinetic reactants
    if (kinetic_reactants)
    {
        for (std::size_t chemical_system_id = 0;
             chemical_system_id < num_chemical_systems;
             ++chemical_system_id)
        {
            kinetic_reactants->setChemicalSystemID(chemical_system_id);
            Rxn_kinetics_map.emplace(chemical_system_id,
                                     *kinetic_reactants->castToBaseClass());
        }
        // explicitly release and delete kinetic_reactants
        kinetic_reactants.reset(nullptr);

        use.Set_kinetics_in(true);
    }

    // rates
    reinitializeRates();

    setConvergenceTolerance();

    configureOutputSettings();

    for (auto const& map_pair : process_id_to_component_name_map)
    {
        auto const transport_process_id = map_pair.first;
        auto const& transport_process_variable = map_pair.second;

        auto master_species =
            master_bsearch(transport_process_variable.c_str());

        _process_id_to_master_map[transport_process_id] = master_species;
    }
}

void PhreeqcKernel::tidyEquilibriumReactants(
    EquilibriumReactants const equilibrium_reactants)
{
    // extract a part of function body from int
    // Phreeqc::tidy_pp_assemblage(void)
    count_elts = 0;
    double coef = 1.0;
    for (auto const& phase_component :
         equilibrium_reactants.getPhaseComponents())
    {
        int phase_id;
        struct phase* phase_component_ptr =
            phase_bsearch(phase_component.first.c_str(), &phase_id, FALSE);
        add_elt_list(phase_component_ptr->next_elt, coef);
    }

    cxxNameDouble nd = elt_list_NameDouble();
    const_cast<cxxPPassemblage*>(equilibrium_reactants.castToBaseClass())
        ->Set_eltList(nd);
}

void PhreeqcKernel::loadDatabase(std::string const& database)
{
    std::ifstream in(database);
    if (!in)
    {
        OGS_FATAL("Unable to open database file '{:s}'.", database);
    }
    assert(phrq_io->get_istream() == nullptr);
    phrq_io->push_istream(&in, false);
    read_database();
}

void PhreeqcKernel::reinitializeRates()
{
    count_rates = _reaction_rates.size();
    rates = (struct rate*)realloc(
        rates, (std::size_t)(count_rates) * sizeof(struct rate));
    int rate_id = 0;
    for (auto const& reaction_rate : _reaction_rates)
    {
        rates[rate_id].name = reaction_rate.kinetic_reactant.data();

        // Commands' strings are freed by Phreeqc dtor.
        rates[rate_id].commands = static_cast<char*>(
            malloc(sizeof(char) * reaction_rate.commands().size() + 1));
        if (rates[rate_id].commands == nullptr)
        {
            OGS_FATAL("Could not allocate memory for rate[{:d}] commands.",
                      rate_id);
        }
        reaction_rate.commands().copy(rates[rate_id].commands,
                                      std::string::npos);
        rates[rate_id].commands[reaction_rate.commands().size()] = '\0';

        rates[rate_id].new_def = 1;
        rates[rate_id].linebase = nullptr;
        rates[rate_id].varbase = nullptr;
        rates[rate_id].loopbase = nullptr;

        ++rate_id;
    };
}

void PhreeqcKernel::doWaterChemistryCalculation(
    std::vector<GlobalVector> const& /*interpolated_process_solutions*/,
    double const dt)
{
    std::vector<GlobalVector*> process_solutions;

    setAqueousSolutions(process_solutions);

    setTimeStepSize(dt);

    execute(process_solutions);
}

void PhreeqcKernel::setAqueousSolutions(
    std::vector<GlobalVector*> const& process_solutions)
{
    // Loop over chemical systems
    std::size_t const num_chemical_systems = process_solutions[0]->size();
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        auto& aqueous_solution = Rxn_solution_map[chemical_system_id];

        auto initial_aqueous_solution =
            getOrCreateInitialAqueousSolution(aqueous_solution);

        auto& components = initial_aqueous_solution->Get_comps();
        // Loop over transport process id map to retrieve component
        // concentrations from process solutions
        for (auto const& map_pair : _process_id_to_master_map)
        {
            auto const transport_process_id = map_pair.first;
            auto const& master_species = map_pair.second;

            auto& transport_process_solution =
                process_solutions[transport_process_id];

            auto& element_name = master_species->elt->name;
            auto const concentration =
                transport_process_solution->get(chemical_system_id);
            if (isHydrogen(element_name))
            {
                // Set pH value by hydrogen concentration.
                double const pH = -std::log10(concentration);
                aqueous_solution.Set_ph(pH);

                {
                    auto hydrogen = components.find("H(1)");
                    if (hydrogen != components.end())
                    {
                        hydrogen->second.Set_input_conc(pH);
                    }
                }
            }
            else
            {
                // Set component concentrations.
                components[element_name].Set_input_conc(concentration);
            }
        }
    }
}

cxxISolution* PhreeqcKernel::getOrCreateInitialAqueousSolution(
    cxxSolution& aqueous_solution)
{
    if (!aqueous_solution.Get_initial_data())
    {
        aqueous_solution.Set_initial_data(_initial_aqueous_solution.get());
        aqueous_solution.Set_new_def(true);
    }

    return aqueous_solution.Get_initial_data();
}

void PhreeqcKernel::setTimeStepSize(double const dt)
{
    // Loop over rxn kinetics map
    for (auto& map_pair : Rxn_kinetics_map)
    {
        auto& kinetics = map_pair.second;
        kinetics.Get_steps().push_back(dt);
    }
}

void PhreeqcKernel::execute(std::vector<GlobalVector*>& process_solutions)
{
    std::size_t const num_chemical_systems = process_solutions[0]->size();
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        Rxn_new_solution.insert(chemical_system_id);
        new_solution = 1;
        use.Set_n_solution_user(chemical_system_id);

        if (!Rxn_pp_assemblage_map.empty())
        {
            Rxn_new_pp_assemblage.insert(chemical_system_id);
            use.Set_pp_assemblage_in(true);
            use.Set_n_pp_assemblage_user(chemical_system_id);
        }

        if (!Rxn_kinetics_map.empty())
        {
            use.Set_kinetics_in(true);
            use.Set_n_kinetics_user(chemical_system_id);
        }

        initial_solutions(false);

        reactions();

        updateNodalProcessSolutions(process_solutions, chemical_system_id);

        reset(chemical_system_id);
    }
}

void PhreeqcKernel::reset(std::size_t const chemical_system_id)
{
    // Clean up
    Rxn_new_solution.clear();

    // Solution
    {
        Rxn_solution_map[chemical_system_id] = *_aqueous_solution;
        Rxn_solution_map[chemical_system_id].Set_n_user_both(
            chemical_system_id);
        Rxn_solution_map[chemical_system_id].Set_pe(-s_eminus->la);
    }

    // Equilibrium reactants
    if (!Rxn_pp_assemblage_map.empty())
    {
        Rxn_new_pp_assemblage.clear();
        for (int j = 0; j < count_unknowns; j++)
        {
            if (x == nullptr || x[j]->type != PP)
                continue;

            auto& phase_components = Rxn_pp_assemblage_map[chemical_system_id]
                                         .Get_pp_assemblage_comps();
            phase_components[x[j]->pp_assemblage_comp_name].Set_moles(
                x[j]->moles);
        }
    }

    // Kinetics
    if (!Rxn_kinetics_map.empty())
    {
        Rxn_kinetics_map[chemical_system_id].Get_steps().clear();
    }
}

void PhreeqcKernel::executeInitialCalculation(
    std::vector<GlobalVector> const& /*interpolated_process_solutions*/)
{
    // TODO (Renchao): This function could be replaced with
    // PhreeqcKernel::doWaterChemistryCalculation(std::vector<GlobalVector*>&
    // process_solutions, double const dt).
    std::vector<GlobalVector*> process_solutions;

    setAqueousSolutions(process_solutions);

    setTimeStepSize(0);

    execute(process_solutions);
}

void PhreeqcKernel::updateNodalProcessSolutions(
    std::vector<GlobalVector*> const& process_solutions,
    std::size_t const node_id)
{
    for (auto const& map_pair : _process_id_to_master_map)
    {
        auto const transport_process_id = map_pair.first;
        auto const& master_species = map_pair.second;

        auto& transport_process_solution =
            process_solutions[transport_process_id];

        auto const& element_name = master_species->elt->name;
        if (isHydrogen(element_name))
        {
            // Update hydrogen concentration by pH value.
            auto const hydrogen_concentration = std::pow(10, s_hplus->la);
            transport_process_solution->set(node_id, hydrogen_concentration);
        }
        else
        {
            // Update solutions of component transport processes.
            auto const concentration =
                master_species->primary
                    ? master_species->total_primary / mass_water_aq_x
                    : master_species->total / mass_water_aq_x;
            transport_process_solution->set(node_id, concentration);
        }
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
