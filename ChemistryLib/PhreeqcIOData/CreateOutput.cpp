/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateOutput.h"

#include <algorithm>
#include <numeric>

#include "AqueousSolution.h"
#include "ChemicalSystem.h"
#include "KineticReactant.h"
#include "UserPunch.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<Output> createOutput(
    ChemicalSystem const& chemical_system,
    std::unique_ptr<UserPunch> const& user_punch,
    bool const use_high_precision,
    std::string const& project_file_name)
{
    auto const& components = chemical_system.aqueous_solution->components;
    auto const& equilibrium_reactants = chemical_system.equilibrium_reactants;
    auto const& kinetic_reactants = chemical_system.kinetic_reactants;
    // Mark which phreeqc output items will be held.
    std::vector<OutputItem> accepted_items{{"pH", ItemType::pH},
                                           {"pe", ItemType::pe}};

    auto accepted_item = [](auto const& item)
    { return OutputItem(item.name, item.item_type); };
    std::transform(components.begin(), components.end(),
                   std::back_inserter(accepted_items), accepted_item);
    std::transform(equilibrium_reactants.begin(), equilibrium_reactants.end(),
                   std::back_inserter(accepted_items), accepted_item);
    for (auto const& kinetic_reactant : kinetic_reactants)
    {
        if (kinetic_reactant.fix_amount)
        {
            continue;
        }
        accepted_items.emplace_back(kinetic_reactant.name,
                                    kinetic_reactant.item_type);
    }

    if (user_punch)
    {
        auto const& secondary_variables = user_punch->secondary_variables;
        accepted_items.reserve(accepted_items.size() +
                               secondary_variables.size());
        std::transform(secondary_variables.begin(), secondary_variables.end(),
                       std::back_inserter(accepted_items), accepted_item);
    }

    // Record ids of which phreeqc output items will be dropped.
    BasicOutputSetups basic_output_setups(project_file_name,
                                          use_high_precision);
    auto const num_dropped_basic_items =
        basic_output_setups.getNumberOfDroppedItems();
    std::vector<int> dropped_item_ids(num_dropped_basic_items);
    std::iota(dropped_item_ids.begin(), dropped_item_ids.end(), 0);

    auto const num_dvalue_equilibrium_reactants = equilibrium_reactants.size();
    auto const num_dvalue_kinetic_reactants = kinetic_reactants.size();
    int const num_dvalue_items =
        num_dvalue_equilibrium_reactants + num_dvalue_kinetic_reactants;

    auto const num_basic_items =
        basic_output_setups.getNumberOfItemsInDisplay();
    auto const num_components = components.size();
    auto dvalue_item_id = num_basic_items + num_components + 1;
    for (int i = 0; i < num_dvalue_items; ++i)
    {
        dropped_item_ids.push_back(dvalue_item_id);
        dvalue_item_id += 2;
    }

    return std::make_unique<Output>(std::move(basic_output_setups),
                                    std::move(accepted_items),
                                    std::move(dropped_item_ids));
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
