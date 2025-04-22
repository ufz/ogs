/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <string>
#include <string_view>
#include <vector>

#pragma once

namespace MeshLib
{
/// Description of the stored integration point data providing additional
/// information for reconstruction and post-processing.
struct IntegrationPointMetaData
{
    /// Converts integration point meta data to a JSON string.
    static std::string toJsonString(
        std::vector<IntegrationPointMetaData> const& meta_data);

    /// Constructs integration point meta data from a JSON encoded string and
    /// looks for the name to be extracted from the parsed integration point
    /// meta data. Fails if no meta data was found for the given name.
    static IntegrationPointMetaData fromJsonString(
        std::string_view const json_string, std::string const& name);

    std::string name;
    int n_components;
    int integration_order;
};

}  // namespace MeshLib
