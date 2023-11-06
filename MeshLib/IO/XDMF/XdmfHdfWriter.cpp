/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "XdmfHdfWriter.h"

#include <algorithm>
#include <functional>
#include <range/v3/algorithm/contains.hpp>

#include "BaseLib/Algorithm.h"
#include "InfoLib/GitInfo.h"
#include "partition.h"
#include "transformData.h"
#include "writeXdmf.h"

using namespace std::literals;

namespace MeshLib::IO
{
struct TransformedMeshData final
{
    std::vector<double> flattened_geometry_values;
    std::vector<int> flattened_topology_values;
};
struct XdmfHdfMesh final
{
    XdmfHdfData geometry;
    XdmfHdfData topology;
    std::vector<XdmfHdfData> attributes;
    std::string name;
    // TransformedMeshData may be large, ensure it is never copied
    std::unique_ptr<TransformedMeshData> transformed_data;
};

template <typename Data>
std::function<bool(Data)> isVariableAttribute(
    std::set<std::string> const& variable_output_names)
{
    if (variable_output_names.empty())
    {
        return [](Data const& data) -> bool
        {
            constexpr std::array constant_output_names{
                "MaterialIDs"sv,
                "topology"sv,
                "geometry"sv,
                "OGS_VERSION"sv,
                MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
                MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
                MeshLib::getBulkIDString(MeshLib::MeshItemType::Edge),
                MeshLib::getBulkIDString(MeshLib::MeshItemType::Face)};
            return !ranges::contains(constant_output_names, data.name);
        };
    }
    return [&variable_output_names](Data const& data) -> bool
    { return variable_output_names.contains(data.name); };
}

XdmfHdfWriter::XdmfHdfWriter(
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
    std::filesystem::path const& filepath, unsigned long long const time_step,
    double const initial_time,
    std::set<std::string> const& variable_output_names,
    bool const use_compression, unsigned int const n_files,
    unsigned int const chunk_size_bytes)
{
    // ogs meshes to vector of Xdmf/HDF meshes (we keep Xdmf and HDF together
    // because XDMF depends on HDF) to meta

    // if no output name is specified, all data will be assumened to be
    // variable over the timesteps. The xdmfhdfwriter is an alternative to
    // other writers, that do not consider the constantness of data Callers
    // of xdmfwriter (e.g. ogs tools) do not provide these information yet
    // and indicate with empty list

    // Transform the data to be written into a format conforming with the rules
    // of xdmf topology and geometry
    auto const transform_ogs_mesh_data_to_xdmf_conforming_data =
        [&n_files, &chunk_size_bytes](auto const& mesh)
    {
        auto flattened_geometry_values = transformToXDMFGeometry(mesh);
        // actually this line is only needed to calculate the offset
        XdmfHdfData const& geometry = transformGeometry(
            mesh, flattened_geometry_values.data(), n_files, chunk_size_bytes);
        auto const flattened_topology_values =
            transformToXDMFTopology(mesh, geometry.hdf.offsets[0]);
        return std::make_unique<TransformedMeshData>(
            TransformedMeshData{std::move(flattened_geometry_values),
                                std::move(flattened_topology_values)});
    };

    // create metadata for transformed data and original ogs mesh data
    auto const transform_to_meta_data =
        [&transform_ogs_mesh_data_to_xdmf_conforming_data, &n_files,
         &chunk_size_bytes](auto const& mesh)
    {
        // important: transformed data must survive and be unique, raw pointer
        // to its memory!
        std::unique_ptr<TransformedMeshData> xdmf_conforming_data =
            transform_ogs_mesh_data_to_xdmf_conforming_data(mesh);
        auto const geometry = transformGeometry(
            mesh, xdmf_conforming_data->flattened_geometry_values.data(),
            n_files, chunk_size_bytes);
        auto const topology =
            transformTopology(xdmf_conforming_data->flattened_topology_values,
                              n_files, chunk_size_bytes);
        auto const attributes =
            transformAttributes(mesh, n_files, chunk_size_bytes);
        return XdmfHdfMesh{std::move(geometry), std::move(topology),
                           std::move(attributes), mesh.get().getName(),
                           std::move(xdmf_conforming_data)};
    };
    auto isVariableHdfAttribute =
        isVariableAttribute<HdfData>(variable_output_names);

    // extract meta data relevant for HDFWriter
    auto const transform_metamesh_to_hdf =
        [&isVariableHdfAttribute](auto const& metamesh)
    {
        // topology and geometry can be treated as any other attribute
        std::vector<HdfData> hdf_data_attributes = {metamesh.geometry.hdf,
                                                    metamesh.topology.hdf};

        hdf_data_attributes.reserve(hdf_data_attributes.size() +
                                    metamesh.attributes.size());
        std::transform(metamesh.attributes.begin(), metamesh.attributes.end(),
                       std::back_inserter(hdf_data_attributes),
                       [](XdmfHdfData att) -> HdfData { return att.hdf; });

        HDFAttributes constant_attributes;
        std::copy_if(hdf_data_attributes.begin(), hdf_data_attributes.end(),
                     back_inserter(constant_attributes),
                     std::not_fn(isVariableHdfAttribute));
        HDFAttributes variable_attributes;
        std::copy_if(hdf_data_attributes.begin(), hdf_data_attributes.end(),
                     back_inserter(variable_attributes),
                     isVariableHdfAttribute);

        return MeshHdfData{
            .constant_attributes = std::move(constant_attributes),
            .variable_attributes = std::move(variable_attributes),
            .name = std::move(metamesh.name)};
    };

    // --------------- XDMF + HDF ---------------------
    std::vector<XdmfHdfMesh> xdmf_hdf_meshes;
    xdmf_hdf_meshes.reserve(meshes.size());
    std::transform(meshes.begin(), meshes.end(),
                   std::back_inserter(xdmf_hdf_meshes), transform_to_meta_data);

    std::vector<MeshHdfData> hdf_meshes;
    hdf_meshes.reserve(xdmf_hdf_meshes.size());
    std::transform(xdmf_hdf_meshes.begin(), xdmf_hdf_meshes.end(),
                   std::back_inserter(hdf_meshes), transform_metamesh_to_hdf);

    // --------------- HDF ---------------------
    std::filesystem::path const hdf_filepath =
        filepath.parent_path() / (filepath.stem().string() + ".h5");

    auto const is_file_manager = isFileManager();
    _hdf_writer = std::make_unique<HdfWriter>(std::move(hdf_meshes), time_step,
                                              hdf_filepath, use_compression,
                                              is_file_manager, n_files);

    // --------------- XDMF ---------------------
    // The light data is only written by just one process
    if (!is_file_manager)
    {
        return;
    }

    auto isVariableXdmfAttribute =
        isVariableAttribute<XdmfData>(variable_output_names);
    // xdmf section
    // extract meta data relevant for XDMFWriter
    auto const transform_metamesh_to_xdmf =
        [&isVariableXdmfAttribute, &filepath, &hdf_filepath,
         &initial_time](XdmfHdfMesh& metamesh)
    {
        std::string const xdmf_name = metamesh.name;
        std::filesystem::path const xdmf_filepath =
            filepath.parent_path() /
            (filepath.stem().string() + "_" + xdmf_name + ".xdmf");

        std::vector<XdmfData> xdmf_attributes;
        std::transform(metamesh.attributes.begin(), metamesh.attributes.end(),
                       std::back_inserter(xdmf_attributes),
                       [](XdmfHdfData const& att) -> XdmfData
                       { return att.xdmf; });

        for (std::size_t i = 0; i < metamesh.attributes.size(); ++i)
        {
            // index 1 time,  index 2 geo, index 3 topology, attributes start at
            // index 4
            xdmf_attributes[i].index = i + 4;
        }

        std::vector<XdmfData> xdmf_variable_attributes;
        std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                     back_inserter(xdmf_variable_attributes),
                     isVariableXdmfAttribute);
        std::vector<XdmfData> xdmf_constant_attributes;
        std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                     back_inserter(xdmf_constant_attributes),
                     std::not_fn(isVariableXdmfAttribute));

        auto const xdmf_writer_fn =
            write_xdmf(metamesh.geometry.xdmf, metamesh.topology.xdmf,
                       xdmf_constant_attributes, xdmf_variable_attributes,
                       hdf_filepath.filename().string(),
                       GitInfoLib::GitInfo::ogs_version, xdmf_name);
        auto xdmf_writer = std::make_unique<XdmfWriter>(xdmf_filepath.string(),
                                                        xdmf_writer_fn);
        xdmf_writer->addTimeStep(initial_time);
        return xdmf_writer;
    };

    std::transform(xdmf_hdf_meshes.begin(), xdmf_hdf_meshes.end(),
                   std::back_inserter(_xdmf_writer),
                   transform_metamesh_to_xdmf);
}

void XdmfHdfWriter::writeStep(double const time)
{
    // ToDo (tm) time_step will be used for simulation continuation (restart)
    _hdf_writer->writeStep(time);
    // The light data is only written by just one process
    if (isFileManager())
    {
        for (auto const& xdmf_writer : _xdmf_writer)
        {
            xdmf_writer->addTimeStep(time);
        }
    }
}

}  // namespace MeshLib::IO
