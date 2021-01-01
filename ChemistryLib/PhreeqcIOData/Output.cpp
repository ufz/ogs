/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "Output.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::ostream& operator<<(std::ostream& os,
                         BasicOutputSetups const& basic_output_setups)
{
    os << "SELECTED_OUTPUT" << "\n";
    os << "-file " << basic_output_setups.output_file << "\n";
    os << "-high_precision " << std::boolalpha
       << basic_output_setups.use_high_precision << "\n";
    os << "-simulation " << std::boolalpha
       << BasicOutputSetups::display_simulation_id << "\n";
    os << "-state " << std::boolalpha << BasicOutputSetups::display_state
       << "\n";
    os << "-distance " << std::boolalpha << BasicOutputSetups::display_distance
       << "\n";
    os << "-time " << std::boolalpha << BasicOutputSetups::display_current_time
       << "\n";
    os << "-step " << std::boolalpha << BasicOutputSetups::display_time_step
       << "\n";

    return os;
}

std::ostream& operator<<(std::ostream& os, Output const& output)
{
    os << output.basic_output_setups;

    auto const component_items =
        output.getOutputItemsByItemType(ItemType::Component);
    os << "-totals";
    for (auto const& component_item : component_items)
    {
        os << " " << component_item.name;
    }
    os << "\n";

    auto const equilibrium_phase_items =
        output.getOutputItemsByItemType(ItemType::EquilibriumReactant);
    if (!equilibrium_phase_items.empty())
    {
        os << "-equilibrium_phases";
        for (auto const& equilibrium_phase_item : equilibrium_phase_items)
        {
            os << " " << equilibrium_phase_item.name;
        }
        os << "\n";
    }

    auto const kinetic_reactant_items =
        output.getOutputItemsByItemType(ItemType::KineticReactant);
    if (!kinetic_reactant_items.empty())
    {
        os << "-kinetic_reactants";
        for (auto const& kinetic_reactant_item : kinetic_reactant_items)
        {
            os << " " << kinetic_reactant_item.name;
        }
        os << "\n";
    }

    return os;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
