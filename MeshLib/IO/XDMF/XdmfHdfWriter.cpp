// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    std::vector<std::size_t> flattened_topology_values;
    ParentDataType parent_data_type;
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

constexpr std::array constant_output_names{
    "MaterialIDs"sv,
    "topology"sv,
    "geometry"sv,
    "OGS_VERSION"sv,
    MeshLib::getBulkIDString(MeshLib::MeshItemType::Node),
    MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell),
    MeshLib::getBulkIDString(MeshLib::MeshItemType::Edge),
    MeshLib::getBulkIDString(MeshLib::MeshItemType::Face),
    MeshLib::globalIDString(MeshLib::MeshItemType::Node),
    MeshLib::globalIDString(MeshLib::MeshItemType::Cell)};

bool isStaticAttribute(std::string const& attribute_name,
                       std::set<std::string> const& static_attribute_names)
{
    return ranges::contains(constant_output_names, attribute_name) ||
           static_attribute_names.contains(attribute_name);
}

template <typename Data>
std::function<bool(Data)> isVariableAttribute(
    std::set<std::string> const& output_variable_names,
    std::set<std::string> const& static_attribute_names)
{
    if (output_variable_names.empty())
    {
        return [&static_attribute_names](Data const& data) -> bool
        { return !isStaticAttribute(data.name, static_attribute_names); };
    }
    return [&output_variable_names,
            &static_attribute_names](Data const& data) -> bool
    {
        if (isStaticAttribute(data.name, static_attribute_names))
        {
            return false;
        }
        return output_variable_names.contains(data.name);
    };
}

XdmfHdfWriter::XdmfHdfWriter(
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> const& meshes,
    std::filesystem::path const& filepath, unsigned long long const time_step,
    double const initial_time,
    std::set<std::string> const& output_variable_names,
    std::set<std::string> const& static_attribute_names,
    bool const use_compression, unsigned int const n_files,
    unsigned int const chunk_size_bytes,
    bool const store_static_data_separately)
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
        auto const [flattened_topology_values, parent_data_type] =
            transformToXDMFTopology(mesh, geometry.hdf.offsets[0]);
        return std::make_unique<TransformedMeshData>(TransformedMeshData{
            std::move(flattened_geometry_values),
            std::move(flattened_topology_values), parent_data_type});
    };

    // create metadata for transformed data and original ogs mesh data
    auto const transform_to_meta_data =
        [&output_variable_names,
         &transform_ogs_mesh_data_to_xdmf_conforming_data, &n_files,
         &chunk_size_bytes](auto const& mesh)
    {
        // important: transformed data must survive and be unique, raw pointer
        // to its memory!
        std::unique_ptr<TransformedMeshData> xdmf_conforming_data =
            transform_ogs_mesh_data_to_xdmf_conforming_data(mesh);
        auto const geometry = transformGeometry(
            mesh, xdmf_conforming_data->flattened_geometry_values.data(),
            n_files, chunk_size_bytes);
        auto const topology = transformTopology(
            xdmf_conforming_data->flattened_topology_values,
            xdmf_conforming_data->parent_data_type, n_files, chunk_size_bytes);
        auto const attributes = transformAttributes(mesh, output_variable_names,
                                                    n_files, chunk_size_bytes);
        return XdmfHdfMesh{std::move(geometry), std::move(topology),
                           std::move(attributes), mesh.get().getName(),
                           std::move(xdmf_conforming_data)};
    };
    auto isVariableHdfAttribute = isVariableAttribute<HdfData>(
        output_variable_names, static_attribute_names);

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
    std::filesystem::path const static_hdf_filepath =
        filepath.parent_path() / (filepath.stem().string() + "_static.h5");

    auto const is_file_manager = isFileManager();

    if (store_static_data_separately)
    {
        // Split hdf meshes into static (constant attributes only) and dynamic
        // (variable attributes only).
        std::vector<MeshHdfData> static_hdf_meshes;
        std::vector<MeshHdfData> dynamic_hdf_meshes;
        for (auto& m : hdf_meshes)  // move data out of m
        {
            static_hdf_meshes.push_back(MeshHdfData{
                .constant_attributes = std::move(m.constant_attributes),
                .variable_attributes = {},
                .name = m.name});
            dynamic_hdf_meshes.push_back(MeshHdfData{
                .constant_attributes = {},
                .variable_attributes = std::move(m.variable_attributes),
                .name = std::move(m.name)});
        }

        _static_hdf_writer = std::make_unique<HdfWriter>(
            std::move(static_hdf_meshes), time_step, initial_time,
            static_hdf_filepath, use_compression, is_file_manager, n_files,
            true);
        _hdf_writer = std::make_unique<HdfWriter>(
            std::move(dynamic_hdf_meshes), time_step, initial_time,
            hdf_filepath, use_compression, is_file_manager, n_files, false);
    }
    else
    {
        // Backward-compatible single-file mode: all data goes to one file with
        // time dimension.
        _hdf_writer = std::make_unique<HdfWriter>(
            std::move(hdf_meshes), time_step, initial_time, hdf_filepath,
            use_compression, is_file_manager, n_files, false);
    }

    // --------------- XDMF ---------------------
    // The light data is only written by just one process
    if (!is_file_manager)
    {
        return;
    }

    auto isVariableXdmfAttribute = isVariableAttribute<XdmfData>(
        output_variable_names, static_attribute_names);
    // xdmf section
    // extract meta data relevant for XDMFWriter
    auto const transform_metamesh_to_xdmf =
        [&isVariableXdmfAttribute, &filepath, &hdf_filepath,
         &static_hdf_filepath, &initial_time,
         store_static_data_separately](XdmfHdfMesh& metamesh)
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

        std::vector<XdmfData> xdmf_variable_attributes;
        std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                     back_inserter(xdmf_variable_attributes),
                     isVariableXdmfAttribute);
        std::vector<XdmfData> xdmf_constant_attributes;
        std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                     back_inserter(xdmf_constant_attributes),
                     std::not_fn(isVariableXdmfAttribute));

        // In backward compat mode, constant data also has a time dimension and
        // shares the dynamic file.
        auto const static_fn = store_static_data_separately
                                   ? static_hdf_filepath.filename().string()
                                   : hdf_filepath.filename().string();
        auto const xdmf_writer_fn =
            write_xdmf(metamesh.geometry.xdmf, metamesh.topology.xdmf,
                       xdmf_constant_attributes, xdmf_variable_attributes,
                       static_fn, hdf_filepath.filename().string(),
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
