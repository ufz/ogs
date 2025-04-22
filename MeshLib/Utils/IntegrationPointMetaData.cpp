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

#include "BaseLib/Error.h"

using nlohmann::json;

namespace MeshLib
{
std::string IntegrationPointMetaData::toJsonString(
    std::vector<MeshLib::IntegrationPointMetaData> const& meta_data)
{
    json json_meta_data;
    json_meta_data["integration_point_arrays"] = json::array();

    for (auto const& md : meta_data)
    {
        json_meta_data["integration_point_arrays"].push_back(
            {{"name", md.name},
             {"number_of_components", md.n_components},
             {"integration_order", md.integration_order}});
    }

    return json_meta_data.dump();
}

IntegrationPointMetaData IntegrationPointMetaData::fromJsonString(
    std::string_view const json_string, std::string const& name)
{
    json const meta_data = json::parse(json_string);

    // Find the current integration point data entry and extract the
    // meta data.
    auto const& ip_meta_data = meta_data["integration_point_arrays"];
    if (auto const it =
            find_if(cbegin(ip_meta_data), cend(ip_meta_data),
                    [&name](auto const& md) { return md["name"] == name; });
        it != cend(ip_meta_data))
    {
        return {name, (*it)["number_of_components"],
                (*it)["integration_order"]};
    }
    OGS_FATAL("No integration point meta data with name '{:s}' found.", name);
}
}  // namespace MeshLib
