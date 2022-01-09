/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LookupTable.h"

#include <unordered_set>

#include "BaseLib/Algorithm.h"

namespace ProcessLib
{
namespace ComponentTransport
{
static std::vector<std::size_t> intersection(
    std::vector<std::size_t> const& vec1, std::vector<std::size_t> const& vec2)
{
    std::unordered_set<std::size_t> set(vec1.begin(), vec1.end());
    std::vector<std::size_t> vec;
    for (auto const a : vec2)
    {
        if (set.contains(a))
        {
            vec.push_back(a);
            set.erase(a);
        }
    }
    return vec;
}

std::pair<double, double> Field::getBoundingSeedPoints(double const value) const
{
    auto it = std::lower_bound(seed_points.cbegin(), seed_points.cend(), value);
    if (it == seed_points.begin())
    {
        WARN("The interpolation point is below the lower bound.");
        auto const nx = std::next(it);
        return std::make_pair(*it, *nx);
    }
    if (it == seed_points.end())
    {
        WARN("The interpolation point is above the upper bound.");
        std::advance(it, -1);
    }

    auto const pv = std::prev(it);
    return std::make_pair(*pv, *it);
}

void LookupTable::lookup(std::vector<GlobalVector*> const& x,
                         std::vector<GlobalVector*> const& x_prev,
                         std::size_t const n_nodes) const
{
    using EntryInput = std::vector<double>;

    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
    {
        std::vector<InterpolationPoint> interpolation_points;
        EntryInput base_entry_input;
        for (auto const& input_field : input_fields)
        {
            // process id and variable id are equilvalent in the case the
            // staggered coupling scheme is adopted.
            auto const process_id = input_field.variable_id;
            auto const& variable_name = input_field.name;
            double input_field_value =
                variable_name.find("_prev") == std::string::npos
                    ? x[process_id]->get(node_id)
                    : x_prev[process_id]->get(node_id);
            input_field_value =
                (std::abs(input_field_value) + input_field_value) / 2;
            auto bounding_seed_points =
                input_field.getBoundingSeedPoints(input_field_value);

            InterpolationPoint interpolation_point{bounding_seed_points,
                                                   input_field_value};
            interpolation_points.push_back(interpolation_point);

            base_entry_input.push_back(bounding_seed_points.first);
        }

        auto const base_entry_id = getTableEntryID(base_entry_input);

        // collect bounding entry ids
        EntryInput bounding_entry_input{base_entry_input};
        std::vector<std::size_t> bounding_entry_ids;
        for (std::size_t i = 0; i < interpolation_points.size(); ++i)
        {
            bounding_entry_input[i] =
                interpolation_points[i].bounding_points.second;
            bounding_entry_ids.push_back(getTableEntryID(bounding_entry_input));
            bounding_entry_input[i] =
                interpolation_points[i].bounding_points.first;
        }

        for (auto const& input_field : input_fields)
        {
            if (input_field.name.find("_prev") != std::string::npos)
            {
                continue;
            }

            auto const output_field_name = input_field.name + "_new";
            auto const base_value =
                tabular_data.at(output_field_name)[base_entry_id];
            auto new_value = base_value;

            // linear interpolation
            for (std::size_t i = 0; i < interpolation_points.size(); ++i)
            {
                auto const interpolation_point_value =
                    tabular_data.at(output_field_name)[bounding_entry_ids[i]];
                auto const slope =
                    (interpolation_point_value - base_value) /
                    (interpolation_points[i].bounding_points.second -
                     interpolation_points[i].bounding_points.first);

                new_value +=
                    slope * (interpolation_points[i].value -
                             interpolation_points[i].bounding_points.first);
            }

            x[input_field.variable_id]->set(node_id, new_value);
        }
    }
}

std::size_t LookupTable::getTableEntryID(
    std::vector<double> const& entry_input) const
{
    std::vector<std::size_t> temp_vec;

    std::vector<std::size_t> intersected_vec =
        /// point_id_groups stores indices where the elements equal to the given
        /// value.
        input_fields[0].point_id_groups[BaseLib::findIndex(
            input_fields[0].seed_points, entry_input[0])];

    /// need to intersect point_id_groups with respect to different input fields
    /// to find out an entry id where all the field data align with the entry
    /// inputs.
    for (std::size_t i = 1; i < input_fields.size(); ++i)
    {
        std::vector<std::size_t> const vec =
            input_fields[i].point_id_groups[BaseLib::findIndex(
                input_fields[i].seed_points, entry_input[i])];

        temp_vec = intersection(intersected_vec, vec);

        std::swap(intersected_vec, temp_vec);
        temp_vec.clear();
    }

    return intersected_vec[0];
}
}  // namespace ComponentTransport
}  // namespace ProcessLib
