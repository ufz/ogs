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
        OGS_FATAL("Could not open file '%s' for writing phreeqc inputs.",
                  _phreeqc_input_file.c_str());

    out << *this;

    if (!out)
        OGS_FATAL("Failed in generating phreeqc input file '%s'.",
                  _phreeqc_input_file.c_str());

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
}  // namespace ChemistryLib
