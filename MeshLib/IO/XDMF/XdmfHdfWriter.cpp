/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "XdmfHdfWriter.h"

#include <algorithm>
#include <functional>

#include "InfoLib/GitInfo.h"
#include "partition.h"
#include "transformData.h"
#include "writeXdmf.h"

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

XdmfHdfWriter::XdmfHdfWriter(
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> meshes,
    std::filesystem::path const& filepath, unsigned long long const time_step,
    double const initial_time,
    std::set<std::string> const& variable_output_names,
    bool const use_compression)
{
    // ogs meshes to vector of Xdmf/HDF meshes (we keep Xdmf and HDF together
    // because XDMF depends on HDF) to meta

    // if no output name is specified, all data will be assumened to be
    // variable over the timesteps. The xdmfhdfwriter is an alternative to
    // other writers, that do not consider the constantness of data Callers
    // of xdmfwriter (e.g. ogs tools) do not provide these information yet
    // and indicate with empty list
    std::function<bool(HdfData)> const is_variable_hdf_attribute =
        [&variable_output_names](
            bool outputnames) -> std::function<bool(HdfData)>
    {
        if (outputnames)
        {
            return [&variable_output_names](HdfData const& data) -> bool
            {
                return std::find(variable_output_names.begin(),
                                 variable_output_names.end(),
                                 data.name) != variable_output_names.end();
            };
        }
        else
        {
            return [](HdfData const&) -> bool { return true; };
        }
    }(!variable_output_names.empty());

    auto is_variable_xdmf_attribute =
        [&variable_output_names](XdmfData const& data) -> bool
    {
        return std::find(variable_output_names.begin(),
                         variable_output_names.end(),
                         data.name) != variable_output_names.end();
    };

    // Transform the data to be written into a format conforming with the rules
    // of xdmf topology and geometry
    auto const transform_ogs_mesh_data_to_xdmf_conforming_data =
        [](auto const& mesh)
    {
        auto flattened_geometry_values = transformToXDMFGeometry(mesh);
        // actually this line is only needed to calculate the offset
        XdmfHdfData const& geometry =
            transformGeometry(mesh, flattened_geometry_values.data());
        auto const flattened_topology_values =
            transformToXDMFTopology(mesh, geometry.hdf.offsets[0]);
        return std::make_unique<TransformedMeshData>(
            TransformedMeshData{std::move(flattened_geometry_values),
                                std::move(flattened_topology_values)});
    };

    // create metadata for transformed data and original ogs mesh data
    auto const transform_to_meta_data =
        [&transform_ogs_mesh_data_to_xdmf_conforming_data](auto const& mesh)
    {
        // important: transformed data must survive and be unique, raw pointer
        // to its memory!
        std::unique_ptr<TransformedMeshData> xdmf_conforming_data =
            transform_ogs_mesh_data_to_xdmf_conforming_data(mesh);
        auto const geometry = transformGeometry(
            mesh, xdmf_conforming_data->flattened_geometry_values.data());
        auto const topology =
            transformTopology(xdmf_conforming_data->flattened_topology_values);
        auto const attributes = transformAttributes(mesh);
        return XdmfHdfMesh{std::move(geometry), std::move(topology),
                           std::move(attributes), mesh.get().getName(),
                           std::move(xdmf_conforming_data)};
    };

    // extract meta data relevant for HDFWriter
    auto const transform_metamesh_to_hdf =
        [&is_variable_hdf_attribute](auto const& metamesh)
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
                     std::not_fn(is_variable_hdf_attribute));
        HDFAttributes variable_attributes;
        std::copy_if(hdf_data_attributes.begin(), hdf_data_attributes.end(),
                     back_inserter(variable_attributes),
                     is_variable_hdf_attribute);

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
                                              is_file_manager);

    // --------------- XDMF ---------------------
    // The light data is only written by just one process
    if (!is_file_manager)
    {
        return;
    }

    // xdmf section
    // extract meta data relevant for XDMFWriter
    auto const transform_metamesh_to_xdmf =
        [&is_variable_xdmf_attribute, &filepath, &hdf_filepath,
         &initial_time](XdmfHdfMesh& metamesh)
    {
        std::string const xdmf_name = metamesh.name;
        std::filesystem::path const xdmf_filepath =
            filepath.parent_path() / (xdmf_name + ".xdmf");

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
                     is_variable_xdmf_attribute);
        std::vector<XdmfData> xdmf_constant_attributes;
        std::copy_if(xdmf_attributes.begin(), xdmf_attributes.end(),
                     back_inserter(xdmf_constant_attributes),
                     std::not_fn(is_variable_xdmf_attribute));

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
