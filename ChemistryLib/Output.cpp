/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "Output.h"

namespace ChemistryLib
{
std::ofstream& operator<<(std::ofstream& out,
                          BasicOutputSetups const& basic_output_setups)
{
    out << "-file " << basic_output_setups.output_file << "\n";
    out << "-high_precision " << std::boolalpha
        << basic_output_setups.use_high_precision << "\n";
    out << "-simulation " << std::boolalpha
        << basic_output_setups.display_simulation_id << "\n";
    out << "-state " << std::boolalpha << basic_output_setups.display_state
        << "\n";
    out << "-distance " << std::boolalpha
        << basic_output_setups.display_distance << "\n";
    out << "-time " << std::boolalpha
        << basic_output_setups.display_current_time << "\n";
    out << "-step " << std::boolalpha << basic_output_setups.display_time_step
        << "\n";

    return out;
}

std::ofstream& operator<<(std::ofstream& out, Output const& output)
{
    out << output.basic_output_setups;

    auto const component_items =
        output.getOutputItemsByItemType(ItemType::Component);
    out << "-totals";
    for (auto const& component_item : component_items)
    {
        out << " " << component_item.name;
    }
    out << "\n";

    auto const equilibrium_phase_items =
        output.getOutputItemsByItemType(ItemType::EquilibriumPhase);
    if (!equilibrium_phase_items.empty())
    {
        out << "-equilibrium_phases";
        for (auto const& equilibrium_phase_item : equilibrium_phase_items)
        {
            out << " " << equilibrium_phase_item.name;
        }
        out << "\n";
    }

    auto const kinetic_reactant_items =
        output.getOutputItemsByItemType(ItemType::KineticReactant);
    if (!kinetic_reactant_items.empty())
    {
        out << "-kinetic_reactants";
        for (auto const& kinetic_reactant_item : kinetic_reactant_items)
        {
            out << " " << kinetic_reactant_item.name;
        }
        out << "\n";
    }

    return out;
}
}  // namespace ChemistryLib
