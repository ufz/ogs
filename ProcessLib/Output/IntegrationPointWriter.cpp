/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IntegrationPointWriter.h"

#include <nlohmann/json.hpp>

#include "MeshLib/Mesh.h"

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

    for (const auto& element_ip_values : ip_values)
    {
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

/// For the given json object and the name extract integration point meta data,
/// or fail if no meta data was found for the given name.
static ProcessLib::IntegrationPointMetaData extractIntegrationPointMetaData(
    json const& meta_data, std::string const& name)
{
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

namespace ProcessLib
{
void addIntegrationPointWriter(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer)
{
    std::vector<IntegrationPointMetaData> meta_data;
    meta_data.reserve(size(integration_point_writer));
    transform(cbegin(integration_point_writer), cend(integration_point_writer),
              back_inserter(meta_data),
              [&](auto const& ip_writer)
              { return addIntegrationPointData(mesh, *ip_writer); });
    if (!meta_data.empty())
    {
        addIntegrationPointMetaData(mesh, meta_data);
    }
}

IntegrationPointMetaData getIntegrationPointMetaData(MeshLib::Mesh const& mesh,
                                                     std::string const& name)
{
    if (!mesh.getProperties().existsPropertyVector<char>(
            "IntegrationPointMetaData"))
    {
        OGS_FATAL(
            "Integration point data '{:s}' is present in the vtk field data "
            "but the required 'IntegrationPointMetaData' array is not "
            "available.",
            name);
    }
    auto const& mesh_property_ip_meta_data =
        *mesh.getProperties().template getPropertyVector<char>(
            "IntegrationPointMetaData");

    if (mesh_property_ip_meta_data.getMeshItemType() !=
        MeshLib::MeshItemType::IntegrationPoint)
    {
        OGS_FATAL("IntegrationPointMetaData array must be field data.");
    }

    // Find the current integration point data entry and extract the
    // meta data.
    auto const ip_meta_data = extractIntegrationPointMetaData(
        json::parse(mesh_property_ip_meta_data.begin(),
                    mesh_property_ip_meta_data.end()),
        name);

    return ip_meta_data;
}
}  // namespace ProcessLib
