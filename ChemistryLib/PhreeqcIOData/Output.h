// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
/**
 * \class BasicOutputSetups
 * \brief Controls which built-in PHREEQC columns appear in the output file.
 *
 * This class defines the SELECTED_OUTPUT / USER_PUNCH header for PHREEQC:
 * which standard columns (pH, pe, solution ID, etc.) are requested,
 * whether high precision is used, and what output file name PHREEQC writes.
 *
 * The helper methods getNumberOfItemsInDisplay() and
 * getNumberOfDroppedItems() report how many of those columns will be kept
 * or ignored when OpenGeoSys later parses PHREEQC output.
 */
class BasicOutputSetups
{
public:
    explicit BasicOutputSetups(std::string const& project_file_name,
                               bool const use_high_precision_);

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

/**
 * \enum ItemType
 * \brief Category used to interpret a PHREEQC output column.
 *
 * Distinguishes pH / pe, total component amounts, equilibrium-phase
 * amounts, kinetic-phase amounts, and user-defined secondary variables.
 */
enum class ItemType
{
    pH,
    pe,
    Component,
    EquilibriumReactant,
    KineticReactant,
    SecondaryVariable
};

/**
 * \struct OutputItem
 * \brief One PHREEQC output column that OpenGeoSys will keep.
 *
 * name is the column label from PHREEQC output.
 * item_type tells PhreeqcIO how to map that column back into the
 * chemical state (e.g. pH, component total, mineral amount).
 */
struct OutputItem
{
    OutputItem(std::string name_, ItemType item_type_)
        : name(std::move(name_)), item_type(item_type_)
    {
    }

    std::string const name;
    ItemType const item_type;
};

/**
 * \struct Output
 * \brief Specification of which PHREEQC output columns are imported into
 * OpenGeoSys.
 *
 * PHREEQC writes a table of values (pH, pe, total component amounts,
 * equilibrium phase amounts, kinetic phase amounts, user-defined quantities,
 * etc.) for every \c chemical_system_id. OpenGeoSys does not use all columns.
 *
 * Output defines the subset that will be read back:
 *  - basic_output_setups describes how the PHREEQC SELECTED_OUTPUT /
 *    USER_PUNCH is configured (file name, precision, which built-in
 *    columns like pH and pe are requested, which are dropped),
 *  - accepted_items lists the columns that OpenGeoSys will parse and
 *    map back into its state (AqueousSolution, reactant inventories,
 *    secondary variables),
 *  - dropped_item_ids records PHREEQC columns that are printed but ignored.
 *
 * The ordering in accepted_items matches the column ordering produced
 * by PHREEQC. After each chemistry step, PhreeqcIO iterates over these
 * items, per \c chemical_system_id, to assign the updated chemical state.
 */
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
