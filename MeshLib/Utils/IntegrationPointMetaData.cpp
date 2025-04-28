/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IntegrationPointMetaData.h"

#include <nlohmann/json.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/view/filter.hpp>

#include "BaseLib/Error.h"

using nlohmann::json;

namespace MeshLib
{
std::string IntegrationPointMetaDataSingleField::toJsonString(
    std::vector<MeshLib::IntegrationPointMetaDataSingleField> const& meta_data)
{
    json json_meta_data;
    json_meta_data["integration_point_arrays"] = json::array();

    for (auto const& md : meta_data)
    {
        json_meta_data["integration_point_arrays"].push_back(
            {{"name", md.field_name},
             {"number_of_components", md.n_components},
             {"integration_order", md.integration_order}});
    }

    return json_meta_data.dump();
}

IntegrationPointMetaDataSingleField
IntegrationPointMetaDataSingleField::fromJsonString(
    std::string_view const json_string, std::string const& field_name)
{
    json const meta_data = json::parse(json_string);
    auto const& ip_meta_data = meta_data["integration_point_arrays"];

    auto ip_meta_data_for_name =
        ip_meta_data |
        ranges::views::filter([&field_name](auto const& md)
                              { return md["name"] == field_name; });

    auto const count = ranges::distance(ip_meta_data_for_name);
    if (count == 0)
    {
        OGS_FATAL("No integration point meta data with name '{:s}' found.",
                  field_name);
    }
    if (count != 1)
    {
        OGS_FATAL(
            "Expected exactly one integration point meta data with name "
            "'{:s}', found {}.",
            field_name, count);
    }

    auto const& md = *ranges::begin(ip_meta_data_for_name);
    return {field_name, md["number_of_components"], md["integration_order"]};
}
}  // namespace MeshLib
