/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <numeric>

#include "AqueousSolution.h"
#include "CreateOutput.h"
#include "EquilibriumPhase.h"
#include "KineticReactant.h"
#include "UserPunch.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<Output> createOutput(
    std::vector<Component> const& components,
    std::vector<EquilibriumPhase> const& equilibrium_phases,
    std::vector<KineticReactant> const& kinetic_reactants,
    std::unique_ptr<UserPunch> const& user_punch,
    std::string const& project_file_name)
{
    // Mark which phreeqc output items will be held.
    std::vector<OutputItem> accepted_items{{"pH", ItemType::pH},
                                           {"pe", ItemType::pe}};

    auto accepted_item = [](auto const& item) {
        return OutputItem(item.name, item.item_type);
    };
    std::transform(components.begin(), components.end(),
                   std::back_inserter(accepted_items), accepted_item);
    std::transform(equilibrium_phases.begin(), equilibrium_phases.end(),
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
        std::transform(secondary_variables.begin(), secondary_variables.end(),
                       std::back_inserter(accepted_items), accepted_item);
    }

    // Record ids of which phreeqc output items will be dropped.
    BasicOutputSetups basic_output_setups(project_file_name);
    auto const num_dropped_basic_items =
        basic_output_setups.getNumberOfDroppedItems();
    std::vector<int> dropped_item_ids(num_dropped_basic_items);
    std::iota(dropped_item_ids.begin(), dropped_item_ids.end(), 0);

    auto const num_dvalue_equilibrium_phases = equilibrium_phases.size();
    auto const num_dvalue_kinetic_reactants = kinetic_reactants.size();
    int const num_dvalue_items =
        num_dvalue_equilibrium_phases + num_dvalue_kinetic_reactants;

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
