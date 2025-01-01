/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/adjacent_remove_if.hpp>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/PropertyVector.h"

namespace MaterialLib
{

std::vector<int> splitMaterialIdString(std::string const& material_id_string)
{
    auto const material_ids_strings =
        BaseLib::splitString(material_id_string, ',');

    std::vector<int> material_ids;
    for (auto& mid_str : material_ids_strings)
    {
        std::size_t num_chars_processed = 0;
        int material_id;
        try
        {
            material_id = std::stoi(mid_str, &num_chars_processed);
        }
        catch (std::invalid_argument&)
        {
            OGS_FATAL(
                "Could not parse material ID from '{}' to a valid integer.",
                mid_str);
        }
        catch (std::out_of_range&)
        {
            OGS_FATAL(
                "Could not parse material ID from '{}'. The integer value of "
                "the given string exceeds the permitted range.",
                mid_str);
        }

        if (num_chars_processed != mid_str.size())
        {
            // Not the whole string has been parsed. Check the rest.
            if (auto const it = std::find_if_not(
                    begin(mid_str) + num_chars_processed, end(mid_str),
                    [](unsigned char const c) { return std::isspace(c); });
                it != end(mid_str))
            {
                OGS_FATAL(
                    "Could not parse material ID from '{}'. Please separate "
                    "multiple material IDs by comma only. Invalid character: "
                    "'{}' at position {}.",
                    mid_str, *it, distance(begin(mid_str), it));
            }
        }

        material_ids.push_back(material_id);
    };

    return material_ids;
}

std::vector<int> parseMaterialIdString(
    std::string const& material_id_string,
    MeshLib::PropertyVector<int> const* const material_ids)
{
    if (material_id_string == "*")
    {
        if (material_ids == nullptr)
        {
            OGS_FATAL(
                "MaterialIDs property is not defined in the mesh but it is "
                "required to parse '*' definition.");
        }

        std::vector<int> material_ids_of_this_medium =
            *material_ids |
            ranges::views::adjacent_remove_if(std::equal_to<>()) |
            ranges::to_vector;
        BaseLib::makeVectorUnique(material_ids_of_this_medium);
        DBUG("Catch all medium definition for material ids {}.",
             fmt::join(material_ids_of_this_medium, ", "));
        return material_ids_of_this_medium;
    }

    // Usual case of ids separated by comma.
    return splitMaterialIdString(material_id_string);
}

}  // namespace MaterialLib
