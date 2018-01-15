/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <nlohmann/json.hpp>

#include "MeshLib/Mesh.h"

#include "IntegrationPointWriter.h"

using nlohmann::json;

/// Adds the integration point data and creates meta data for it.
///
/// Returns meta data for the written integration point data.
static ProcessLib::IntegrationPointMetaData addIntegrationPointData(
    MeshLib::Mesh& mesh, ProcessLib::IntegrationPointWriter const& writer)
{
    auto const& ip_values = writer.values(/*t, x, dof_table*/);
    assert(ip_values.size() == mesh.getNumberOfElements());

    // create field data and fill it with nodal values, and an offsets cell
    // array indicating where the cell's integration point data starts.
    auto& field_data = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, writer.name(), MeshLib::MeshItemType::IntegrationPoint,
        writer.numberOfComponents());
    field_data.clear();

    for (std::size_t e = 0; e < ip_values.size(); ++e)
    {
        auto const& element_ip_values = ip_values[e];
        std::copy(element_ip_values.begin(), element_ip_values.end(),
                  std::back_inserter(field_data));
    }

    return {writer.name(), writer.numberOfComponents(),
            writer.integrationOrder()};
}

/// Adds integration point meta data as char mesh property encoded in JSON
/// format, which is then stored as VTK's field data.
static void addIntegrationPointMetaData(
    MeshLib::Mesh& mesh,
    std::vector<ProcessLib::IntegrationPointMetaData> const& meta_data)
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

    // Store the field data.
    std::string const json_string = json_meta_data.dump();
    auto& dictionary = *MeshLib::getOrCreateMeshProperty<char>(
        mesh, "IntegrationPointMetaData",
        MeshLib::MeshItemType::IntegrationPoint, 1);
    dictionary.clear();
    std::copy(json_string.begin(), json_string.end(),
              std::back_inserter(dictionary));
}

namespace ProcessLib
{
void addIntegrationPointWriter(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer)
{
    std::vector<IntegrationPointMetaData> meta_data;
    for (auto const& ip_writer : integration_point_writer)
    {
        meta_data.push_back(addIntegrationPointData(mesh, *ip_writer));
    }
    if (!meta_data.empty())
    {
        addIntegrationPointMetaData(mesh, meta_data);
    }
}
}  // namespace ProcessLib
