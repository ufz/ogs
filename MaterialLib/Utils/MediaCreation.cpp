// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <spdlog/fmt/ranges.h>

#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/unique.hpp>

#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/PropertyVector.h"

namespace MaterialLib
{

/// Checks that a string part contains no whitespace.
/// Throws OGS_FATAL if any whitespace is found.
void checkForWhitespaces(std::string_view part)
{
    if (std::ranges::any_of(
            part,
            [](char c) { return std::isspace(static_cast<unsigned char>(c)); }))
    {
        OGS_FATAL(
            "Whitespace is not allowed in ranges. Use 'start:end' without "
            "spaces around the colon.");
    }
}

/// Creates a range of integers from start to end (inclusive).
/// Throws OGS_FATAL if end < start.
auto expandRange(int start, int end)
{
    if (end < start)
    {
        OGS_FATAL(
            "Invalid range '{}:{}'. The end must be greater than or equal "
            "to the start.",
            start, end);
    }
    return ranges::views::iota(start, end + 1);
}

std::vector<int> splitMaterialIdString(std::string const& material_id_string)
{
    auto const material_ids_strings =
        BaseLib::splitString(material_id_string, ',');

    // Pre-allocate with estimated capacity (simplified heuristic)
    std::vector<int> material_ids;
    material_ids.reserve(material_ids_strings.size());

    for (std::string mid_str : material_ids_strings)
    {
        // Trim leading and trailing whitespace
        BaseLib::trim(mid_str);

        auto const parts =
            BaseLib::splitString(mid_str, ':') | ranges::to_vector;
        if (parts.size() == 2)
        {
            checkForWhitespaces(parts[0]);
            auto const start_id = BaseLib::parseInteger(parts[0]);
            if (!start_id)
            {
                OGS_FATAL("Could not parse material ID: {}", start_id.error());
            }
            checkForWhitespaces(parts[1]);
            auto const end_id = BaseLib::parseInteger(parts[1]);
            if (!end_id)
            {
                OGS_FATAL("Could not parse material ID: {}", end_id.error());
            }
            ranges::copy(expandRange(*start_id, *end_id),
                         std::back_inserter(material_ids));
        }
        else if (parts.size() == 1)
        {
            auto const material_id = BaseLib::parseInteger(mid_str);
            if (!material_id)
            {
                OGS_FATAL("Could not parse material ID: {}",
                          material_id.error());
            }
            material_ids.push_back(*material_id);
        }
        else
        {
            OGS_FATAL(
                "Could not parse material ID from '{}'. Invalid range format. "
                "Use 'start:end' for ranges or a single integer.",
                mid_str);
        }
    }

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

        std::vector<int> const material_ids_of_this_medium =
            *material_ids | ranges::views::unique | ranges::to_vector |
            ranges::actions::sort | ranges::actions::unique |
            ranges::to<std::vector>;
        DBUG("Catch all medium definition for material ids {}.",
             fmt::join(material_ids_of_this_medium, ", "));
        return material_ids_of_this_medium;
    }

    // Usual case of ids or ranges separated by comma.
    return splitMaterialIdString(material_id_string);
}

}  // namespace MaterialLib
