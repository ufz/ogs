/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/cxxKinetics.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
PhreeqcKernel::PhreeqcKernel(
    std::size_t const num_chemical_systems,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    std::string const& database,
    PhreeqcKernelData::AqueousSolution& aqueous_solution,
    std::unique_ptr<Kinetics>
        kinetic_reactants,
    std::vector<ReactionRate>&& reaction_rates)
    : Phreeqc(),
      _initial_aqueous_solution(aqueous_solution.getInitialAqueousSolution()),
      _reaction_rates(std::move(reaction_rates))
{
    do_initialize();

    // load database
    std::ifstream in(database);
    if (!in)
    {
        OGS_FATAL("Unable to open database file '%s'.", database.c_str());
    }
    assert(phrq_io->get_istream() == nullptr);
    phrq_io->push_istream(&in, false);
    read_database();

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

    // kinetics
    if (kinetic_reactants->isKineticReactantDefined())
    {
        for (std::size_t chemical_system_id = 0;
             chemical_system_id < num_chemical_systems;
             ++chemical_system_id)
        {
            kinetic_reactants->setChemicalSystemID(chemical_system_id);
            Rxn_kinetics_map.emplace(chemical_system_id,
                                     *kinetic_reactants->castToBaseClass());
        }
        use.Set_kinetics_in(true);
    }

    reinitializeRates();

    setConvergenceTolerance();

    configureOutputSettings();

    for (auto const& [transport_process_id, transport_process_variable] :
         process_id_to_component_name_map)
    {
        auto master_species =
            master_bsearch(transport_process_variable.c_str());

        _process_id_to_master_map[transport_process_id] = master_species;
    }
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
            OGS_FATAL("Could not allocate memory for rate[%d] commands.",
                      rate_id);
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
    std::vector<GlobalVector*>& process_solutions, double const dt)
{
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
        for (auto const& [transport_process_id, master_species] :
             _process_id_to_master_map)
        {
            auto& transport_process_solution =
                process_solutions[transport_process_id];

            auto& element_name = master_species->elt->name;
            if (isHydrogen(element_name))
            {
                // Set pH value by hydrogen concentration.
                double const pH = -std::log10(
                    transport_process_solution->get(chemical_system_id));
                aqueous_solution.Set_ph(pH);
            }
            else
            {
                // Set component concentrations.
                auto const concentration =
                    transport_process_solution->get(chemical_system_id);
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
        aqueous_solution.Set_initial_data(&_initial_aqueous_solution);
        aqueous_solution.Set_new_def(true);
    }

    return aqueous_solution.Get_initial_data();
}

void PhreeqcKernel::setTimeStepSize(double const dt)
{
    // Loop over rxn kinetics map
    for (auto& [chemical_system_id, kinetics] : Rxn_kinetics_map)
    {
        (void)chemical_system_id;
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

        if (!Rxn_kinetics_map.empty())
        {
            use.Set_kinetics_in(true);
            use.Set_n_kinetics_user(chemical_system_id);
        }

        initial_solutions(false);

        reactions();

        updateNodalProcessSolutions(process_solutions, chemical_system_id);

        // Clean up
        Rxn_new_solution.clear();
        Rxn_solution_map[chemical_system_id].Get_totals().clear();

        if (!Rxn_kinetics_map.empty())
        {
            Rxn_kinetics_map[chemical_system_id].Get_steps().clear();
        }
    }
}

void PhreeqcKernel::updateNodalProcessSolutions(
    std::vector<GlobalVector*> const& process_solutions,
    std::size_t const node_id)
{
    for (auto const& [transport_process_id, master_species] :
         _process_id_to_master_map)
    {
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
                master_species->total_primary / mass_water_aq_x;
            transport_process_solution->set(node_id, concentration);
        }
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
