/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "Output.h"
#include "PhreeqcIO.h"
#include "ThirdParty/iphreeqc/src/src/IPhreeqc.h"

namespace ChemistryLib
{
void PhreeqcIO::doWaterChemistryCalculation(
    std::vector<GlobalVector*>& process_solutions, double const dt)
{
    setAqueousSolutionsOrUpdateProcessSolutions(
        process_solutions, Status::SettingAqueousSolutions);
    setTimeStep(dt);

    writeInputsToFile();

    execute();

    readOutputsFromFile();

    setAqueousSolutionsOrUpdateProcessSolutions(
        process_solutions, Status::UpdatingProcessSolutions);
}

void PhreeqcIO::setAqueousSolutionsOrUpdateProcessSolutions(
    std::vector<GlobalVector*> const& process_solutions, Status const status)
{
    std::size_t const num_chemical_systems = _aqueous_solutions.size();
    // Loop over chemical systems
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        // Get chemical compostion of solution in a particular chemical system
        auto& aqueous_solution = _aqueous_solutions[chemical_system_id];
        auto& components = aqueous_solution.components;
        // Loop over transport process id map to retrieve component
        // concentrations from process solutions or to update process solutions
        // after chemical calculation by Phreeqc
        for (auto const& process_id_to_component_name_map_element :
             _process_id_to_component_name_map)
        {
            auto const& transport_process_id =
                process_id_to_component_name_map_element.first;
            auto const& transport_process_variable =
                process_id_to_component_name_map_element.second;

            auto& transport_process_solution =
                process_solutions[transport_process_id];

            auto component =
                std::find_if(components.begin(), components.end(),
                             [&transport_process_variable](Component const& c) {
                                 return c.name == transport_process_variable;
                             });

            if (component != components.end())
            {
                switch (status)
                {
                    case Status::SettingAqueousSolutions:
                        // Set component concentrations.
                        component->amount =
                            transport_process_solution->get(chemical_system_id);
                        break;
                    case Status::UpdatingProcessSolutions:
                        // Update solutions of component transport processes.
                        transport_process_solution->set(chemical_system_id,
                                                        component->amount);
                        break;
                }
            }

            if (transport_process_variable == "H")
            {
                switch (status)
                {
                    case Status::SettingAqueousSolutions:
                    {
                        // Set pH value by hydrogen concentration.
                        aqueous_solution.pH =
                            -std::log10(transport_process_solution->get(
                                chemical_system_id));
                        break;
                    }
                    case Status::UpdatingProcessSolutions:
                    {
                        // Update hydrogen concentration by pH value.
                        auto hydrogen_concentration =
                            std::pow(10, -aqueous_solution.pH);
                        transport_process_solution->set(chemical_system_id,
                                                        hydrogen_concentration);
                        break;
                    }
                }
            }
        }
    }
}

void PhreeqcIO::writeInputsToFile()
{
    DBUG("Writing phreeqc inputs into file '%s'.", _phreeqc_input_file.c_str());
    std::ofstream out(_phreeqc_input_file, std::ofstream::out);

    if (!out)
    {
        OGS_FATAL("Could not open file '%s' for writing phreeqc inputs.",
                  _phreeqc_input_file.c_str());
    }

    out << *this;

    if (!out)
    {
        OGS_FATAL("Failed in generating phreeqc input file '%s'.",
                  _phreeqc_input_file.c_str());
    }

    out.close();
}

std::ofstream& operator<<(std::ofstream& out, PhreeqcIO const& phreeqc_io)
{
    out << "SELECTED_OUTPUT" << "\n";
    out << *phreeqc_io._output << "\n";

    auto const& reaction_rates = phreeqc_io._reaction_rates;
    if (!reaction_rates.empty())
    {
        out << "RATES" << "\n";
        out << reaction_rates << "\n";
    }

    std::size_t const num_chemical_systems =
        phreeqc_io._aqueous_solutions.size();
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        auto const& aqueous_solution =
            phreeqc_io._aqueous_solutions[chemical_system_id];
        out << "SOLUTION " << chemical_system_id + 1 << "\n";
        out << aqueous_solution << "\n";

        auto const& equilibrium_phases =
            phreeqc_io._equilibrium_phases[chemical_system_id];
        if (!equilibrium_phases.empty())
        {
            out << "EQUILIBRIUM_PHASES " << chemical_system_id + 1 << "\n";
            out << equilibrium_phases << "\n";
        }

        auto const& kinetic_reactants =
            phreeqc_io._kinetic_reactants[chemical_system_id];
        if (!kinetic_reactants.empty())
        {
            out << "KINETICS " << chemical_system_id + 1 << "\n";
            out << kinetic_reactants;
            out << "-steps " << phreeqc_io._dt << "\n" << "\n";
        }

        out << "END" << "\n" << "\n";
    }

    return out;
}

void PhreeqcIO::execute()
{
    auto database_loaded = [&](int instance_id) {
        return LoadDatabase(instance_id, _database.c_str()) == 0;
    };

    INFO("Phreeqc: Executing chemical calculation.");
    // initialize phreeqc configurations
    auto const instance_id = CreateIPhreeqc();

    // load a specific database in the working directory
    if (!database_loaded(instance_id))
    {
        OGS_FATAL(
            "Failed in loading the specified thermodynamic database file: %s.",
            _database.c_str());
    }

    SetSelectedOutputFileOn(instance_id, 1);

    if (RunFile(instance_id, _phreeqc_input_file.c_str()) > 0)
    {
        OutputErrorString(instance_id);
        OGS_FATAL(
            "Failed in performing speciation calculation with the generated "
            "phreeqc input file '%s'.",
            _phreeqc_input_file.c_str());
    }
}

void PhreeqcIO::readOutputsFromFile()
{
    auto const& basic_output_setups = _output->basic_output_setups;
    auto const& phreeqc_result_file = basic_output_setups.output_file;
    DBUG("Reading phreeqc results from file '%s'.",
         phreeqc_result_file.c_str());
    std::ifstream in(phreeqc_result_file);

    if (!in)
    {
        OGS_FATAL("Could not open phreeqc result file '%s'.",
                  phreeqc_result_file.c_str());
    }

    in >> *this;

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc result file '%s'",
                  phreeqc_result_file.c_str());
    }

    in.close();
}

std::ifstream& operator>>(std::ifstream& in, PhreeqcIO& phreeqc_io)
{
    // Skip the headline
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string line;
    auto const& output = *phreeqc_io._output;
    auto const& dropped_item_ids = output.dropped_item_ids;
    std::size_t const num_chemical_systems =
        phreeqc_io._aqueous_solutions.size();
    for (std::size_t chemical_system_id = 0;
         chemical_system_id < num_chemical_systems;
         ++chemical_system_id)
    {
        // Skip equilibrium calculation result of initial solution
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Get calculation result of the solution after the reaction
        if (!std::getline(in, line))
        {
            OGS_FATAL(
                "Error when reading calculation result of Solution %u after "
                "the reaction.",
                chemical_system_id);
        }

        std::vector<double> accepted_items;
        std::vector<std::string> items;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::algorithm::split(items, line, boost::is_any_of("\t "),
                                boost::token_compress_on);
        for (int item_id = 0; item_id < static_cast<int>(items.size());
             ++item_id)
        {
            if (std::find(dropped_item_ids.begin(), dropped_item_ids.end(),
                          item_id) == dropped_item_ids.end())
            {
                accepted_items.push_back(std::stod(items[item_id]));
            }
        }
        assert(accepted_items.size() == output.accepted_items.size());

        auto& aqueous_solution =
            phreeqc_io._aqueous_solutions[chemical_system_id];
        auto& components = aqueous_solution.components;
        auto& equilibrium_phases =
            phreeqc_io._equilibrium_phases[chemical_system_id];
        auto& kinetic_reactants =
            phreeqc_io._kinetic_reactants[chemical_system_id];
        for (int item_id = 0; item_id < static_cast<int>(accepted_items.size());
             ++item_id)
        {
            auto const& accepted_item = output.accepted_items[item_id];
            auto const& item_name = accepted_item.name;
            switch (accepted_item.item_type)
            {
                case ItemType::pH:
                {
                    // Update pH value
                    aqueous_solution.pH = accepted_items[item_id];
                    break;
                }
                case ItemType::pe:
                {
                    // Update pe value
                    aqueous_solution.pe = accepted_items[item_id];
                    break;
                }
                case ItemType::Component:
                {
                    // Update component concentrations
                    auto& component = BaseLib::findElementOrError(
                        components.begin(), components.end(),
                        [&item_name](Component const& component) {
                            return component.name == item_name;
                        },
                        "Could not find component '" + item_name + "'.");
                    component.amount = accepted_items[item_id];
                    break;
                }
                case ItemType::EquilibriumPhase:
                {
                    // Update amounts of equilibrium phases
                    auto& equilibrium_phase = BaseLib::findElementOrError(
                        equilibrium_phases.begin(), equilibrium_phases.end(),
                        [&item_name](
                            EquilibriumPhase const& equilibrium_phase) {
                            return equilibrium_phase.name == item_name;
                        },
                        "Could not find equilibrium phase '" + item_name +
                            "'.");
                    equilibrium_phase.amount = accepted_items[item_id];
                    break;
                }
                case ItemType::KineticReactant:
                {
                    // Update amounts of kinetic reactants
                    auto& kinetic_reactant = BaseLib::findElementOrError(
                        kinetic_reactants.begin(), kinetic_reactants.end(),
                        [&item_name](KineticReactant const& kinetic_reactant) {
                            return kinetic_reactant.name == item_name;
                        },
                        "Could not find kinetic reactant '" + item_name + "'.");
                    kinetic_reactant.amount = accepted_items[item_id];
                    break;
                }
            }
        }
    }

    return in;
}
}  // namespace ChemistryLib
