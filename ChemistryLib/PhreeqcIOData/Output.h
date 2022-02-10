/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <iosfwd>
#include <iterator>
#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
class BasicOutputSetups
{
public:
    explicit BasicOutputSetups(std::string const& project_file_name,
                               bool const use_high_precision_)
        : output_file(project_file_name + "_phreeqc.out"),
          use_high_precision(use_high_precision_)
    {
    }

    static int getNumberOfItemsInDisplay()
    {
        return display_simulation_id + display_state + display_solution_id +
               display_distance + display_current_time + display_time_step +
               display_pH + display_pe;
    }

    static int getNumberOfDroppedItems()
    {
        return display_simulation_id + display_state + display_solution_id +
               display_distance + display_current_time + display_time_step;
    }

    friend std::ostream& operator<<(
        std::ostream& os, BasicOutputSetups const& basic_output_setups);

    std::string const output_file;

private:
    static const bool display_simulation_id = false;
    static const bool display_state = true;
    static const bool display_solution_id = true;
    static const bool display_distance = false;
    static const bool display_current_time = false;
    static const bool display_time_step = false;
    static const bool display_pH = true;
    static const bool display_pe = true;
    bool const use_high_precision;
};

enum class ItemType
{
    pH,
    pe,
    Component,
    EquilibriumReactant,
    KineticReactant,
    SecondaryVariable
};

struct OutputItem
{
    OutputItem(std::string name_, ItemType item_type_)
        : name(std::move(name_)), item_type(item_type_)
    {
    }

    std::string const name;
    ItemType const item_type;
};

struct Output
{
    Output(BasicOutputSetups&& basic_output_setups_,
           std::vector<OutputItem>&& accepted_items_,
           std::vector<int>&& dropped_item_ids_)
        : basic_output_setups(std::move(basic_output_setups_)),
          accepted_items(std::move(accepted_items_)),
          dropped_item_ids(std::move(dropped_item_ids_))
    {
    }

    std::vector<OutputItem> getOutputItemsByItemType(ItemType item_type) const
    {
        std::vector<OutputItem> matching_items;
        std::copy_if(accepted_items.cbegin(),
                     accepted_items.cend(),
                     std::back_inserter(matching_items),
                     [&item_type](OutputItem const& accepted_item)
                     { return accepted_item.item_type == item_type; });
        return matching_items;
    }

    friend std::ostream& operator<<(std::ostream& os, Output const& output);

    BasicOutputSetups const basic_output_setups;
    std::vector<OutputItem> const accepted_items;
    std::vector<int> const dropped_item_ids;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
