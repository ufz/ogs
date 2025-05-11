/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "IntegrationPointWriter.h"

#include <range/v3/view/join.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"

/// Adds the integration point data and creates meta data for it.
///
/// Returns meta data for the written integration point data.
static MeshLib::IntegrationPointMetaDataSingleField addIntegrationPointData(
    MeshLib::Mesh& mesh, MeshLib::IntegrationPointWriter const& writer)
{
    auto const& ip_values = writer.values(/*t, x, dof_table*/);
    assert(ip_values.size() == mesh.getNumberOfElements());

    // create field data and fill it with nodal values, and an offsets cell
    // array indicating where the cell's integration point data starts.
    auto& field_data = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, writer.name(), MeshLib::MeshItemType::IntegrationPoint,
        writer.numberOfComponents());
    field_data.clear();

    ranges::copy(ip_values | ranges::views::join,
                 std::back_inserter(field_data));

    return {writer.name(), writer.numberOfComponents(),
            writer.integrationOrder()};
}

/// Adds integration point meta data as char mesh property encoded in JSON
/// format, which is then stored as VTK's field data.
static void addIntegrationPointMetaDataSingleField(
    MeshLib::Mesh& mesh, MeshLib::IntegrationPointMetaData const& ip_meta_data)
{
    // Store the field data.
    std::string const json_string = ip_meta_data.toJsonString();
    auto& dictionary = *MeshLib::getOrCreateMeshProperty<char>(
        mesh, "IntegrationPointMetaData",
        MeshLib::MeshItemType::IntegrationPoint, 1);
    dictionary.clear();
    std::copy(json_string.begin(), json_string.end(),
              std::back_inserter(dictionary));
}

namespace MeshLib
{
void addIntegrationPointDataToMesh(
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer)
{
    auto meta_data = IntegrationPointMetaData{
        integration_point_writer |
        ranges::views::transform(
            [&](auto const& ip_writer)
            { return addIntegrationPointData(mesh, *ip_writer); })};

    if (!meta_data.empty())
    {
        addIntegrationPointMetaDataSingleField(mesh, meta_data);
    }
}

std::optional<IntegrationPointMetaData> getIntegrationPointMetaData(
    MeshLib::Properties const& properties)
{
    if (!properties.existsPropertyVector<char>("IntegrationPointMetaData"))
    {
        return std::nullopt;
    }
    auto const& mesh_property_ip_meta_data =
        *properties.template getPropertyVector<char>(
            "IntegrationPointMetaData");

    if (mesh_property_ip_meta_data.getMeshItemType() !=
        MeshLib::MeshItemType::IntegrationPoint)
    {
        OGS_FATAL("IntegrationPointMetaData array must be field data.");
    }

    // Find the current integration point data entry and extract the
    // meta data.
    return IntegrationPointMetaData{std::string_view{
        mesh_property_ip_meta_data.data(), mesh_property_ip_meta_data.size()}};
}
}  // namespace MeshLib
