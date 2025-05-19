/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IntegrationPointMetaData.h"

#include <spdlog/fmt/ranges.h>

#include <nlohmann/json.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <unordered_set>

#include "BaseLib/Error.h"

using nlohmann::json;

namespace MeshLib
{

IntegrationPointMetaData::IntegrationPointMetaData(
    std::string_view const json_string)
{
    json const meta_data = json::parse(json_string);

    fields_ = meta_data["integration_point_arrays"] |
              ranges::views::transform(
                  [](auto const& md)
                  {
                      return IntegrationPointMetaDataSingleField{
                          md["name"], md["number_of_components"],
                          md["integration_order"]};
                  }) |
              ranges::to<std::vector>;

    checkFieldNamesAreUnique();
}

IntegrationPointMetaDataSingleField const& IntegrationPointMetaData::operator[](
    std::string const& field_name) const
{
    auto it = ranges::find(
        fields_, field_name, &IntegrationPointMetaDataSingleField::field_name);

    if (it == fields_.end())
    {
        OGS_FATAL("No integration point meta data with name '{:s}' found.",
                  field_name);
    }
    return *it;
}

std::string IntegrationPointMetaData::toJsonString() const
{
    json json_meta_data;
    json_meta_data["integration_point_arrays"] = json::array();

    for (auto const& field : fields_)
    {
        json_meta_data["integration_point_arrays"].push_back(
            {{"name", field.field_name},
             {"number_of_components", field.n_components},
             {"integration_order", field.integration_order}});
    }

    return json_meta_data.dump();
}

void IntegrationPointMetaData::checkFieldNamesAreUnique() const
{
    std::unordered_set<std::string> seen;

    auto const duplicates =
        fields_ |
        ranges::views::transform(
            &IntegrationPointMetaDataSingleField::field_name) |
        ranges::views::filter([&seen](std::string const& name)
                              { return !seen.insert(name).second; }) |
        ranges::to<std::vector>;

    if (!duplicates.empty())
    {
        OGS_FATAL("Duplicate integration point meta data names found: {:s}.",
                  fmt::join(duplicates, ", "));
    }
}

IntegrationPointMetaDataSingleField getIntegrationPointMetaDataSingleField(
    std::optional<IntegrationPointMetaData> const& ip_meta_data,
    std::string const& field_name)
{
    if (!ip_meta_data)
    {
        OGS_FATAL(
            "The required 'IntegrationPointMetaData' array is not available in "
            "the vtk field data but is needed to evaluate the integration "
            "point property '{:s}'.",
            field_name);
    }
    return (*ip_meta_data)[field_name];
}
}  // namespace MeshLib
