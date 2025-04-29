/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <concepts>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/traits.hpp>
#if defined(__cpp_lib_containers_ranges) && \
    __cpp_lib_containers_ranges >= 202202L
#include <ranges>
#else
#include <range/v3/to_container.hpp>
#endif
#include <string>
#include <string_view>
#include <vector>

#pragma once

namespace MeshLib
{
/// Single field integration point meta data providing additional information
/// for reconstruction and post-processing.
struct IntegrationPointMetaDataSingleField
{
    std::string field_name;
    int n_components;
    int integration_order;

    auto operator<=>(IntegrationPointMetaDataSingleField const&) const =
        default;
};

/// Description of the stored integration point meta data for all fields.
class IntegrationPointMetaData
{
public:
    /// Constructs integration point meta data from a JSON encoded string.
    explicit IntegrationPointMetaData(std::string_view const json_string);

    explicit IntegrationPointMetaData(ranges::input_range auto&& range)
        requires(std::convertible_to<ranges::range_reference_t<decltype(range)>,
                                     IntegrationPointMetaDataSingleField>)
        :
#if defined(__cpp_lib_containers_ranges) && \
    __cpp_lib_containers_ranges >= 202202L
          fields_(std::from_range, std::forward<decltype(range)>(range))
#else
          fields_(ranges::to<std::vector<IntegrationPointMetaDataSingleField>>(
              std::forward<decltype(range)>(range)))
#endif
    {
        checkFieldNamesAreUnique();
    }
    /// Converts integration point meta data to a JSON string.
    std::string toJsonString() const;

    constexpr bool empty() const { return fields_.empty(); }

    IntegrationPointMetaDataSingleField const& operator[](
        std::string const& field_name) const;

    auto operator<=>(IntegrationPointMetaData const&) const = default;

private:
    void checkFieldNamesAreUnique() const;

private:
    std::vector<IntegrationPointMetaDataSingleField> fields_;
};

/// Returns integration point meta data for single field or an error if
/// an empty optional ip_meta_data is passed.
IntegrationPointMetaDataSingleField getIntegrationPointMetaDataSingleField(
    std::optional<IntegrationPointMetaData> const& ip_meta_data,
    std::string const& field_name);
}  // namespace MeshLib
